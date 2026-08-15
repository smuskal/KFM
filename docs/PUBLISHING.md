# Publishing a model bundle

Two commands. Nothing else is a supported way to publish.

```bash
python3 scripts/publish-bundle.py --bundle-dir ./LigASeqLigB_v2_potency          # dry run
python3 scripts/publish-bundle.py --bundle-dir ./LigASeqLigB_v2_potency --apply  # publish
python3 scripts/verify-published-bundles.py --weights                            # before announcing
```

Dry run is the default. Nothing reaches the bucket without `--apply`.

## Why the script exists

Publishing used to be manual. On 2026-08-13 the potency bundle's `MANIFEST.json`
was overwritten with a copy of the bundle's file list, because `manifest.json`
and `MANIFEST.json` are one file on a case-insensitive filesystem. Every
`install.sh` run failed for two days. Nobody noticed, because nothing checked.
The same collision had already destroyed the selectivity manifest once before.

There was no script to fix. That was the actual problem.

## What the script refuses to do

| It stops when | Because |
|---|---|
| Two filenames differ only in case | That pair is one file locally and two objects in the bucket |
| A `manifest.json` sits in the bundle directory | It is generated. Writing it there is what destroys `MANIFEST.json` |
| `MANIFEST.json` holds a file list | That is the exact 2026-08-13 damage |
| `feature_order` is neither `[ ligand` nor `[ sequence` | The tools identify a model from it |
| `model_sha256` does not match the `.joblib` present | The manifest describes a different model |
| Anything fails to verify after upload | A publish that does not verify is a failed publish |

## How it publishes

Each file is uploaded individually to an explicit destination. No rsync, no
recursive copy, no reliance on a local name matching a remote one. The file list
is generated into a temporary directory under a name that cannot collide with
anything in the bundle, then uploaded under its real name.

Data files go first and the file list last, so the list never advertises a file
that has not landed.

Weights get `max-age=604800, immutable` because their filenames carry a date and
are never rewritten. Everything else gets `no-cache`, so a repair is visible at
once instead of sitting behind an hour of cache looking like it failed.

After uploading, every file is re-downloaded and checked against the list. Pass
`--verify-weights` to include the model files.

## If a publish fails verification

Object versioning is on for `gs://kfm-models-public`, so the previous generation
of every object is recoverable:

```bash
gcloud storage ls -a gs://kfm-models-public/<bundle>/
gcloud storage cp gs://kfm-models-public/<bundle>/<file>#<generation> \
                  gs://kfm-models-public/<bundle>/<file>
```

Versioning was off until 2026-08-15, which is why the original potency manifest
had to be reconstructed rather than restored.

## Still worth doing

Rename the file list to `files.json` so the collision cannot exist at all.
`kfm/bundles.py` requests `manifest.json`, so it needs a release that serves both
names before the old one is retired.
