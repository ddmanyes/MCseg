import { copyFile, mkdir, readFile, writeFile } from 'node:fs/promises'
import { createHash } from 'node:crypto'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..')
const artifact = path.resolve(process.argv[2] || '')
if (!/\.(dmg|exe)$/.test(artifact)) throw new Error('Pass the newly built .dmg or .exe path')
const manifestPath = path.join(root, 'frontend/src-tauri/build-resources/app/BUILD_MANIFEST.json')
const manifest = JSON.parse(await readFile(manifestPath, 'utf8'))
if (manifest.source_dirty) throw new Error('Refusing release metadata for dirty sources')
const name = path.basename(artifact)
if (!name.includes(`_${manifest.desktop_version}_`)) throw new Error('Artifact version does not match manifest')
const bytes = await readFile(artifact)
const digest = createHash('sha256').update(bytes).digest('hex')
const output = path.join(root, 'release')
await mkdir(output, {recursive:true})
const metadataName = name + '.build.json'
await copyFile(artifact, path.join(output, name))
await writeFile(path.join(output, metadataName), JSON.stringify({
  ...manifest, artifact: {name, bytes: bytes.length, sha256: digest},
  validation: 'Build and payload provenance only; does not certify clean-machine installation or scientific results',
}, null, 2) + '\n')
await writeFile(path.join(output, name + '.sha256'), `${digest}  ${name}\n`)
console.log(`Finalized ${name}; source ${manifest.source_commit}; SHA256 ${digest}`)
