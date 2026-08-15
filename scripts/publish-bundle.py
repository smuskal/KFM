#!/usr/bin/env python3
"""Publish a model bundle to the public bucket. The only supported way to do it.

    python3 scripts/publish-bundle.py --bundle-dir ./LigASeqLigB_v2_potency
    python3 scripts/publish-bundle.py --bundle-dir ./LigASeqLigB_v2_potency --apply

Dry run by default. It prints exactly what it would upload and exits without
touching the bucket. Nothing is published until you pass --apply.

WHY THIS EXISTS

Publishing used to be manual. On 2026-08-13 the potency bundle's MANIFEST.json
was overwritten with a copy of the bundle's file list, because `manifest.json`
and `MANIFEST.json` are one file on a case-insensitive filesystem and something
staged both from one directory. Every install failed for two days and nothing
caught it. The same collision had already destroyed the selectivity manifest
once before. There was no script to fix, which was the real problem.

WHAT MAKES THAT IMPOSSIBLE HERE

  1. The bundle directory is REJECTED if two filenames differ only in case.
     That is the fault itself, caught before anything is read.
  2. MANIFEST.json is REJECTED if it looks like a file list. That is the exact
     shape of the 2026-08-13 damage, so the specific accident cannot be
     published even if it somehow reaches the staging directory.
  3. Every file is uploaded individually to an explicit destination. No rsync,
     no recursive copy, no reliance on the local name matching the remote one.
  4. The file list is written to a temporary path named after the bundle, never
     alongside a MANIFEST.json, so generating it cannot clobber anything.
  5. After uploading, every file is re-downloaded and checked against the list.
     A publish that does not verify is a failed publish, and it says so.

Requires the gcloud CLI and write access to the bucket.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import urllib.request

BUCKET = "kfm-models-public"
PUBLIC = f"https://storage.googleapis.com/{BUCKET}"
LISTING = "manifest.json"          # the file list the CLI asks for first
BIG = 50_000_000

LICENCE = ("RESEARCH AND EVALUATION USE ONLY. Free for academic and non-profit "
           "research and teaching, and for internal evaluation of up to 30 days. "
           "Commercial use, service provision, redistribution, reverse engineering, "
           "and training any model or building any database from the weights or "
           "their outputs each require a separate written licence.")
LICENCE_URL = "https://kinasefoundationmodel.com/LICENSE-MODELS.txt"
COMMERCIAL = "https://eidogen-sertanty.com/kinasekbmarvin.php"

# Weights carry a date in the filename and are never rewritten, so they may be
# cached hard. Everything else can change between releases and must not be, or a
# repair sits behind the cache for an hour looking like it failed.
CACHE_WEIGHTS = "public, max-age=604800, immutable"
CACHE_OTHER = "no-cache, max-age=0"


def die(msg):
    print(f"\nREFUSING TO PUBLISH: {msg}", file=sys.stderr)
    sys.exit(1)


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def inspect(bundle_dir):
    """Read the bundle and refuse anything that could republish the 13 Aug fault."""
    if not os.path.isdir(bundle_dir):
        die(f"{bundle_dir} is not a directory")
    names = [f for f in sorted(os.listdir(bundle_dir))
             if os.path.isfile(os.path.join(bundle_dir, f)) and not f.startswith(".")]

    # 1. The fault itself. Two names differing only in case are ONE file on a
    # case-insensitive filesystem and TWO objects in the bucket. On Linux both
    # can exist, and MANIFEST.json plus manifest.json is precisely the pair that
    # caused the 2026-08-13 damage. Compared on the real on-disk names, so it
    # fires wherever the pair genuinely exists and not on macOS, where the
    # bundle's own MANIFEST.json is the only file present.
    folded = {}
    for n in names:
        folded.setdefault(n.lower(), []).append(n)
    clashes = {k: v for k, v in folded.items() if len(v) > 1}
    if clashes:
        die("two filenames differ only in case, so they are one file on a "
            f"case-insensitive filesystem: {clashes}. Rename before publishing.")

    # A generated file list has no business in the source directory. Publishing
    # one would ship a stale listing, and on macOS writing it here is what
    # destroys MANIFEST.json. This script always generates it into a temporary
    # directory instead.
    if LISTING in names:
        die(f"{bundle_dir} contains a {LISTING}. That file is generated, and "
            "writing it beside MANIFEST.json is what broke the potency bundle. "
            "Delete it; this script builds the list itself.")

    weights = [n for n in names if n.endswith(".joblib")]
    if len(weights) != 1:
        die(f"a bundle must hold exactly one .joblib; found {weights}")

    mpath = os.path.join(bundle_dir, "MANIFEST.json")
    if not os.path.isfile(mpath):
        die("no MANIFEST.json in the bundle")
    try:
        manifest = json.load(open(mpath))
    except Exception as e:                                   # noqa: BLE001
        die(f"MANIFEST.json is not valid JSON: {e}")

    # 2. The specific damage: a file list published as the model manifest.
    if isinstance(manifest.get("files"), list) and "bundle" in manifest:
        die("MANIFEST.json holds a FILE LIST, not a model manifest. This is "
            "exactly the 2026-08-13 fault. Restore the real manifest first.")
    fo = str(manifest.get("feature_order", ""))
    if not (fo.startswith("[ ligand") or fo.startswith("[ sequence")):
        die(f"MANIFEST.json feature_order is {fo!r}. The tools identify a model "
            "from it, and neither shape matches.")

    declared = manifest.get("model_sha256")
    actual = sha256_of(os.path.join(bundle_dir, weights[0]))
    if declared and declared != actual:
        die(f"MANIFEST.json model_sha256 does not match {weights[0]}.\n"
            f"  manifest {declared}\n  actual   {actual}")
    return names, weights[0], manifest


def build_listing(bundle, bundle_dir, names, weights):
    files = []
    for n in sorted(names):
        p = os.path.join(bundle_dir, n)
        files.append({"name": n, "bytes": os.path.getsize(p), "sha256": sha256_of(p)})
    return {"bundle": bundle, "files": files, "licence": LICENCE,
            "licence_url": LICENCE_URL, "commercial": COMMERCIAL,
            "model_file": weights}


def upload(local, dest, cache):
    subprocess.run(["gcloud", "storage", "cp", local, dest,
                    f"--cache-control={cache}", "--quiet"], check=True)


def verify(bundle, listing, weights_too):
    print("\n--- verifying what is actually served ---")
    bad = 0
    for f in listing["files"] + [{"name": LISTING}]:
        n = f["name"]
        if f.get("bytes", 0) > BIG and not weights_too:
            print(f"  skipped  {n} (pass --verify-weights)")
            continue
        h, size = hashlib.sha256(), 0
        try:
            with urllib.request.urlopen(f"{PUBLIC}/{bundle}/{n}", timeout=1800) as r:
                for chunk in iter(lambda: r.read(1 << 20), b""):
                    h.update(chunk)
                    size += len(chunk)
        except Exception as e:                               # noqa: BLE001
            print(f"  FAILED   {n}: {e}")
            bad += 1
            continue
        want = f.get("sha256")
        if want is None:                                     # the listing itself
            print(f"  ok       {n} (served, {size} bytes)")
        elif h.hexdigest() == want and size == f["bytes"]:
            print(f"  ok       {n}")
        else:
            bad += 1
            print(f"  MISMATCH {n}\n             want {f['bytes']} {want}"
                  f"\n             got  {size} {h.hexdigest()}")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle-dir", required=True)
    ap.add_argument("--name", help="bucket prefix. Default: the directory name")
    ap.add_argument("--apply", action="store_true", help="actually upload")
    ap.add_argument("--verify-weights", action="store_true",
                    help="re-download the weights during verification too")
    a = ap.parse_args()

    bundle_dir = os.path.abspath(a.bundle_dir)
    bundle = a.name or os.path.basename(bundle_dir.rstrip("/"))
    names, weights, _ = inspect(bundle_dir)
    listing = build_listing(bundle, bundle_dir, names, weights)

    print(f"bundle      : {bundle}")
    print(f"source      : {bundle_dir}")
    print(f"weights     : {weights}")
    print(f"files       : {len(names)} plus the generated {LISTING}")
    for f in listing["files"]:
        cache = "immutable" if f["name"] == weights else "no-cache"
        print(f"  {f['bytes']:>12,}  {f['sha256'][:12]}  {f['name']}  [{cache}]")

    if not a.apply:
        print(f"\nDry run. Nothing was uploaded. Re-run with --apply to publish to "
              f"gs://{BUCKET}/{bundle}/")
        return 0

    print(f"\n--- uploading to gs://{BUCKET}/{bundle}/ ---")
    # Files first, listing last. A client reads the listing first, so publishing
    # it last means it never advertises a file that has not landed yet.
    for f in listing["files"]:
        n = f["name"]
        cache = CACHE_WEIGHTS if n == weights else CACHE_OTHER
        print(f"  {n}")
        upload(os.path.join(bundle_dir, n), f"gs://{BUCKET}/{bundle}/{n}", cache)

    # The listing is written under a name that cannot collide with anything in
    # the bundle, in its own temporary directory, and only then uploaded under
    # its real name. This is the step that broke before.
    with tempfile.TemporaryDirectory() as tmp:
        staged = os.path.join(tmp, f"{bundle}.listing.json")
        with open(staged, "w") as fh:
            json.dump(listing, fh, indent=2)
            fh.write("\n")
        print(f"  {LISTING}")
        upload(staged, f"gs://{BUCKET}/{bundle}/{LISTING}", CACHE_OTHER)

    bad = verify(bundle, listing, a.verify_weights)
    if bad:
        print(f"\n{bad} file(s) do not match after upload. THE PUBLISH FAILED.")
        print("Versioning is on for this bucket, so the previous generation can "
              "be restored:\n  gcloud storage ls -a gs://"
              f"{BUCKET}/{bundle}/")
        return 1
    print("\npublished and verified.")
    if not a.verify_weights:
        print("Weights were not re-downloaded. Before announcing a release, run:")
        print("  python3 scripts/verify-published-bundles.py --weights")
    return 0


if __name__ == "__main__":
    sys.exit(main())
