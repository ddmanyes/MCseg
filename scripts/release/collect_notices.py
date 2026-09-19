"""Collect verbatim installed dependency notices for a specific build environment.

Run after npm ci and cargo metadata --locked --format-version 1 --filter-platform
<TARGET>. Does not claim a legal review of all transitive binary dependencies.
"""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cargo-metadata', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    packages = []
    lock = json.loads((root / 'frontend/package-lock.json').read_text())
    for relative, info in lock['packages'].items():
        if not relative:
            continue
        directory = root / 'frontend' / relative
        if not (directory / 'package.json').exists():
            continue  # optional packages for other platforms are not installed
        meta = json.loads((directory / 'package.json').read_text())
        packages.append(('npm', meta.get('name', relative), meta.get('version', ''), meta.get('license', 'not declared'), directory))
    for info in json.loads(args.cargo_metadata.read_text())['packages']:
        if info['name'] == 'mcseg':
            continue
        packages.append(('Rust', info['name'], info['version'], info.get('license') or 'not declared', Path(info['manifest_path']).parent))
    parts = ['# Third-party notices\n',
             'Generated from the installed npm dependency tree and locked Rust metadata used for the build. '
             'Includes build-time dependencies as well as runtime components. '
             'Verbatim license/notice files are included when present at package root or in a licenses directory. '
             'Regenerate on each target before release; this inventory is not a complete legal audit. '
             'Python packages and model weights are downloaded at first launch, retain their own licensing, '
             'and are not redistributed in this installer.\n',
             '## uv 0.12.5 sidecar\n\nSource: https://github.com/astral-sh/uv/tree/0.12.5 '
             '(MIT OR Apache-2.0). The following are the upstream project licenses; '
             'this inventory does not enumerate every dependency statically linked into uv.\n']
    for p in sorted((root / 'licenses/uv').iterdir()):
        parts.append(f'### {p.name}\n\n```text\n{p.read_text()}\n```\n')
    missing = []
    for ecosystem, name, version, license_id, directory in sorted(packages, key=lambda x: x[:3]):
        candidates = []
        for p in directory.iterdir():
            if p.is_file() and p.name.lower().startswith(('license', 'licence', 'copying', 'notice', 'copyright')):
                candidates.append(p)
            elif p.is_dir() and p.name.lower() in ['licenses', 'licences']:
                candidates.extend(f for f in p.rglob('*') if f.is_file())
        parts.append(f'## {ecosystem}: {name} {version}\n\nDeclared license: {license_id}\n')
        if not candidates:
            missing.append(f'{ecosystem}: {name} {version}')
            parts.append('No standalone license text found in this installed package. Consult its upstream source.\n')
        for p in sorted(candidates):
            parts.append(f'### {p.relative_to(directory)}\n\n```text\n{p.read_text(errors="replace")}\n```\n')
    (root / 'THIRD_PARTY_NOTICES.md').write_text('\n'.join(parts))
    print(f'Collected {len(packages)} dependency entries; {len(missing)} without standalone license text.')
    print('\n'.join(missing))


if __name__ == '__main__':
    main()
