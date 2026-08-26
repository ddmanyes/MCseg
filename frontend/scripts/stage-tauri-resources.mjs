import { cp, copyFile, mkdir, rm, stat } from 'node:fs/promises'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const frontendDir = path.resolve(scriptDir, '..')
const repositoryDir = path.resolve(frontendDir, '..')
const srcTauriDir = path.resolve(frontendDir, 'src-tauri')
const stagingDir = path.resolve(srcTauriDir, 'build-resources', 'app')

if (!stagingDir.startsWith(`${srcTauriDir}${path.sep}`)) {
  throw new Error(`Refusing to clean staging directory outside src-tauri: ${stagingDir}`)
}

const excludedDirectoryNames = new Set(['__pycache__', '.pytest_cache', 'tests'])
const shouldCopyRuntimeFile = (source) => {
  const name = path.basename(source)
  if (excludedDirectoryNames.has(name)) return false
  if (name === '.DS_Store' || name.startsWith('._')) return false
  if (name.endsWith('.pyc') || name.endsWith('.pyo')) return false
  return true
}

const copyRequiredDirectory = async (source, destination, filter) => {
  const sourceStats = await stat(source).catch(() => null)
  if (!sourceStats?.isDirectory()) {
    throw new Error(`Required resource directory is missing: ${source}`)
  }
  await cp(source, destination, { recursive: true, filter })
}

await rm(stagingDir, { recursive: true, force: true })
await mkdir(stagingDir, { recursive: true })

await copyRequiredDirectory(
  path.join(repositoryDir, 'backend'),
  path.join(stagingDir, 'backend'),
  shouldCopyRuntimeFile,
)
await copyRequiredDirectory(
  path.join(repositoryDir, 'config'),
  path.join(stagingDir, 'config'),
  (source) => shouldCopyRuntimeFile(source) && path.basename(source) !== 'state.json',
)
await copyRequiredDirectory(
  path.join(frontendDir, 'dist'),
  path.join(stagingDir, 'frontend', 'dist'),
  shouldCopyRuntimeFile,
)

for (const fileName of ['pyproject.toml', 'uv.lock']) {
  await copyFile(path.join(repositoryDir, fileName), path.join(stagingDir, fileName))
}

console.log(`[tauri:prepare] Runtime resources staged at ${stagingDir}`)
