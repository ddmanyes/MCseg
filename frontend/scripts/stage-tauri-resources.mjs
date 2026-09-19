import { copyFile, mkdir, readdir, readFile, rm, lstat, writeFile } from 'node:fs/promises'
import { execFileSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const hash = data => createHash('sha256').update(data).digest('hex')
const excluded = new Set(['tests', '__pycache__', '.pytest_cache', 'state.json'])
const blocked = name => excluded.has(name) || name.startsWith('.') || /\.(pyc|pyo|map)$/.test(name)

export async function stageResources(repositoryDir) {
  const frontendDir = path.join(repositoryDir, 'frontend')
  const stagingDir = path.join(frontendDir, 'src-tauri', 'build-resources', 'app')
  const git = (...args) => execFileSync('git', args, { cwd: repositoryDir, encoding: 'utf8' }).trim()
  const sourceCommit = git('rev-parse', 'HEAD')
  const dirty = Boolean(git('status', '--porcelain', '--untracked-files=no'))
  if (dirty && process.env.MCSEG_ALLOW_DIRTY_BUILD !== '1') {
    throw new Error('Commit tracked changes before a release build (development only: MCSEG_ALLOW_DIRTY_BUILD=1).')
  }
  const tracked = git('ls-files', '-z').split('\0').filter(Boolean)
  const manifestFiles = []
  await rm(stagingDir, { recursive: true, force: true })
  await mkdir(stagingDir, { recursive: true })
  const copy = async (source, relative) => {
    if (!(await lstat(source)).isFile()) throw new Error(`Not a regular file: ${relative}`)
    const data = await readFile(source)
    const target = path.join(stagingDir, relative)
    await mkdir(path.dirname(target), { recursive: true })
    await copyFile(source, target)
    manifestFiles.push({ path: relative.split(path.sep).join('/'), sha256: hash(data), bytes: data.length })
  }
  for (const relative of tracked) {
    if (relative.split('/').some(blocked)) continue
    if ((relative.startsWith('backend/') && relative.endsWith('.py')) ||
        (relative.startsWith('config/') && relative.endsWith('.yaml'))) {
      await copy(path.join(repositoryDir, relative), relative)
    }
  }
  for (const name of ['pyproject.toml', 'uv.lock', 'LICENSE', 'THIRD_PARTY_NOTICES.md']) {
    if (!tracked.includes(name)) throw new Error(`Required tracked resource missing: ${name}`)
    await copy(path.join(repositoryDir, name), name)
  }
  const pipeline = await readFile(path.join(stagingDir, 'config', 'pipeline.yaml'), 'utf8')
  if (/\/Users\/|\/Volumes\/|[A-Za-z]:[\\/]/.test(pipeline) || !/^rois:\s*\[\]\s*$/m.test(pipeline)) {
    throw new Error('Distributable pipeline.yaml must have no local absolute paths or saved ROIs.')
  }
  const copyDist = async (relative = '') => {
    for (const entry of await readdir(path.join(frontendDir, 'dist', relative), { withFileTypes: true })) {
      if (blocked(entry.name)) continue
      const child = path.join(relative, entry.name)
      if (entry.isSymbolicLink()) throw new Error(`Symlink in frontend build: ${child}`)
      if (entry.isDirectory()) await copyDist(child)
      else if (/\.(html|js|css|svg|png|jpg|webp|ico|woff2?)$/.test(entry.name)) {
        await copy(path.join(frontendDir, 'dist', child), path.join('frontend', 'dist', child))
      }
    }
  }
  await copyDist()
  if (!manifestFiles.some(f => f.path === 'frontend/dist/index.html')) throw new Error('Missing built frontend index.html')
  const app = JSON.parse(await readFile(path.join(frontendDir, 'src-tauri', 'tauri.conf.json'), 'utf8'))
  const sidecarName = process.platform === 'win32'
    ? 'uv-x86_64-pc-windows-msvc.exe'
    : `uv-${process.arch === 'arm64' ? 'aarch64' : 'x86_64'}-apple-darwin`
  const sidecarBytes = await readFile(path.join(frontendDir, 'src-tauri', 'binaries', sidecarName))
  const manifest = {
    schema_version: 1, application: 'MCseg', desktop_version: app.version,
    source_commit: sourceCommit, source_dirty: dirty,
    build_platform: process.platform, build_architecture: process.arch,
    uv_sidecar: { name: sidecarName, sha256: hash(sidecarBytes), expected_version: '0.12.5' },
    runtime_dependencies: 'uv.lock; installed on first launch',
    files: manifestFiles.sort((a, b) => a.path.localeCompare(b.path)),
  }
  await writeFile(path.join(stagingDir, 'BUILD_MANIFEST.json'), JSON.stringify(manifest, null, 2) + '\n')
  return manifest
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const repositoryDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..')
  const manifest = await stageResources(repositoryDir)
  console.log(`[tauri:prepare] ${manifest.files.length} files; source ${manifest.source_commit}; dirty=${manifest.source_dirty}`)
}
