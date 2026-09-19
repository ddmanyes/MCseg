import { test } from 'node:test'
import assert from 'node:assert/strict'
import { mkdtemp, mkdir, writeFile, readFile, rm, symlink } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { execFileSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import { stageResources } from '../stage-tauri-resources.mjs'

async function fixture(t) {
  const root = await mkdtemp(path.join(tmpdir(), 'mcseg-stage-test-'))
  t.after(() => rm(root, { recursive: true, force: true }))
  const put = async (name, data) => {
    await mkdir(path.dirname(path.join(root, name)), { recursive: true })
    await writeFile(path.join(root, name), data)
  }
  const git = (...args) => execFileSync('git', args, { cwd: root, stdio: 'pipe' })
  for (const [name, data] of Object.entries({
    'backend/main.py': 'print("runtime")', 'backend/tests/test_private.py': 'private fixture',
    'backend/__pycache__/secret.pyc': 'cache', 'config/state.json': '{"private":"path"}',
    'config/pipeline.yaml': 'paths:\n  data_root: ""\nrois: []\n',
    'frontend/dist/index.html': '<html>MCseg</html>', 'frontend/dist/.env': 'secret',
    'frontend/src-tauri/tauri.conf.json': '{"version":"0.2.1"}',
    'pyproject.toml': '[project]\nname="mcseg"', 'uv.lock': '# lock',
    'LICENSE': 'MIT fixture', 'THIRD_PARTY_NOTICES.md': 'Notices fixture',
  })) await put(name, data)
  const sidecarName = process.platform === 'win32' ? 'uv-x86_64-pc-windows-msvc.exe'
    : `uv-${process.arch === 'arm64' ? 'aarch64' : 'x86_64'}-apple-darwin`
  await put(`frontend/src-tauri/binaries/${sidecarName}`, 'fixture sidecar')
  git('init', '-q'); git('add', '.'); git('-c','user.name=Test','-c','user.email=test@example.invalid','commit','-qm','fixture')
  await put('backend/untracked.py', 'must not ship')
  return { root, put, git }
}

test('ships tracked runtime only; excludes saved state, tests, caches and untracked files; hashes payload', async t => {
  const {root} = await fixture(t)
  const m = await stageResources(root)
  assert.equal(m.source_dirty, false)
  const paths = m.files.map(f => f.path)
  assert(paths.includes('backend/main.py'))
  assert(paths.includes('LICENSE'))
  assert(paths.includes('THIRD_PARTY_NOTICES.md'))
  assert(!paths.some(p => /state.json|test_private|secret|untracked|\.env/.test(p)))
  for (const f of m.files) {
    const bytes = await readFile(path.join(root, 'frontend/src-tauri/build-resources/app', f.path))
    assert.equal(createHash('sha256').update(bytes).digest('hex'), f.sha256)
  }
})
test('rejects dirty release sources', async t => {
  const {root, put} = await fixture(t)
  await put('backend/main.py', 'changed')
  await assert.rejects(stageResources(root), /Commit tracked changes/)
})
test('rejects tracked local paths and saved ROI defaults', async t => {
  const {root, put, git} = await fixture(t)
  await put('config/pipeline.yaml', 'paths:\n  data_root: /Users/example/data\nrois: []\n')
  git('add','.');git('-c','user.name=Test','-c','user.email=test@example.invalid','commit','-qm','bad path')
  await assert.rejects(stageResources(root), /no local absolute paths/)
})
test('rejects frontend symlinks instead of copying external files', async t => {
  const {root} = await fixture(t)
  await symlink(path.join(root, 'LICENSE'), path.join(root, 'frontend/dist/escape.js'))
  await assert.rejects(stageResources(root), /Symlink/)
})
