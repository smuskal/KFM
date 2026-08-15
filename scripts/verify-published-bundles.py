#!/usr/bin/env python3
"""Check every published bundle against its own file list.

    python3 scripts/verify-published-bundles.py            # all bundles, skip the weights
    python3 scripts/verify-published-bundles.py --weights   # include them, about 4.8 GB
    python3 scripts/verify-published-bundles.py --bundle LigASeqLigB_v2_potency

Needs nothing but network access. No credentials, no SDK, no clone of the
bundles. Anyone can run it, which is the point: a user who suspects a bad
download can check the server rather than guess.

Exits 0 when everything matches and 1 when anything does not.

WHY THIS EXISTS
On 2026-08-13 the potency bundle's MANIFEST.json was overwritten with a copy of
the file list, because `manifest.json` and `MANIFEST.json` are one name on a
case-insensitive filesystem and something staged both from one directory. Every
new install failed from that moment until 2026-08-15. Nothing noticed, because
nothing was checking. The same collision had already destroyed the selectivity
manifest once before.

The check itself takes seconds without --weights. Run it after every publish,
and the answer is known immediately instead of arriving as a user's bug report.
"""
import argparse
import hashlib
import sys
import urllib.request

BASE = "https://kinasefoundationmodel.com/api/model"
BUNDLES = ["LigASeqLigB_v2_potency", "SeqALigSeqB_v2_selectivity", "kfm_v1_potency"]
BIG = 50_000_000


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def digest(url):
    h = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(url, timeout=1800) as r:
        for chunk in iter(lambda: r.read(1 << 20), b""):
            h.update(chunk)
            n += len(chunk)
    return h.hexdigest(), n


def check(bundle, weights):
    import json
    print(f"\n=== {bundle} ===")
    try:
        listing = json.loads(fetch(f"{BASE}/{bundle}/manifest.json"))
    except Exception as e:                                   # noqa: BLE001
        print(f"  CANNOT READ THE FILE LIST: {e}")
        return 1
    bad = 0
    for f in listing["files"]:
        name, want, want_bytes = f["name"], f.get("sha256"), f.get("bytes")
        if want_bytes and want_bytes > BIG and not weights:
            print(f"  skipped  {name}  ({want_bytes/2**20:.0f} MB, pass --weights)")
            continue
        try:
            got, n = digest(f"{BASE}/{bundle}/{name}")
        except Exception as e:                               # noqa: BLE001
            print(f"  FAILED   {name}: {e}")
            bad += 1
            continue
        if got == want and (want_bytes is None or n == want_bytes):
            print(f"  ok       {name}")
        else:
            bad += 1
            print(f"  MISMATCH {name}")
            print(f"             want {want_bytes} bytes  {want}")
            print(f"             got  {n} bytes  {got}")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weights", action="store_true",
                    help="also hash the model files, about 4.8 GB of download")
    ap.add_argument("--bundle", action="append", choices=BUNDLES,
                    help="check one bundle; repeatable. Default: all of them")
    a = ap.parse_args()
    bad = sum(check(b, a.weights) for b in (a.bundle or BUNDLES))
    print()
    if bad:
        print(f"{bad} file(s) do not match. The published bundles are NOT intact.")
        print("Do not assume a bad network link: check the bucket before republishing.")
        return 1
    print("every checked file matches its published checksum.")
    if not a.weights:
        print("The model weights were skipped. Re-run with --weights before a release.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
