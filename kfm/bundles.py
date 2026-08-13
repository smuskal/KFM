"""Finding, downloading and loading a model bundle.

A "bundle" is the directory the model ships as: the forest, the sequence-vector
cache, the index files, and — the part that matters — **the model's own
`predict.py`**. Scoring always goes through that file. This package never
reimplements the maths, because the two models do genuinely different things
with their inputs (the potency forest averages both ligand orders to compensate
for being only approximately antisymmetric; the selectivity forest does not,
because its training augmentation already made it symmetric). A reimplementation
that got that wrong would disagree with every published number while looking
perfectly reasonable.

WHERE THE WEIGHTS COME FROM. They are NOT in the git repository — they are
several gigabytes, GitHub rejects anything over 100 MB, and they are under a
different licence from this code (see LICENSE-MODELS.txt). They are fetched from
kinasefoundationmodel.com on first use and cached.

WHERE THEY ARE PUT. **`./kfm-models`, in the directory you run the download
from**, overridable with `KFM_HOME`. They land where you are working, in plain
sight, not in a hidden cache under your home directory: several gigabytes that
someone cannot find is a support question, and a user who does not know where
their models went cannot delete them either.

Because a project directory is the natural home for them, a later command run
from a SUBdirectory still finds them: resolution walks up from the current
directory looking for an existing `kfm-models`, the way `git` finds `.git`. Only
when nothing is found does it fall back to `./kfm-models` for a fresh download.

One thing to avoid: do not put `KFM_HOME` inside Dropbox, iCloud or Drive. A 3 GB
forest in a synced tree gets uploaded, re-downloaded on every other machine, and
counted against the quota, which is nobody's intention.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import urllib.request

# Where the site serves bundles from. Each download is logged there, which is
# how we know the models are being used and by roughly whom.
BASE_URL = os.environ.get("KFM_DOWNLOAD_BASE",
                          "https://kinasefoundationmodel.com/api/model")

MODELS = {
    # v1 answers a DIFFERENT question from either v2 model: it predicts a pIC50
    # value for one compound against one kinase, rather than comparing two
    # things. Only the validated arm is published -- the frontier arm was fitted
    # on every measurement with nothing held back, so no honest accuracy exists
    # for it, and a downloadable model with no measurable accuracy is not
    # something to hand out.
    "v1": {
        "bundle": "kfm_v1_potency",
        "entry": "kfm_predict.py",
        "weights": "kfm_v1_pointwise_20260809.joblib",
        "ram_gb": 6.4,
        "blurb": "Predicted pIC50 for a compound against one kinase (version 1)",
    },
    "potency": {
        "bundle": "LigASeqLigB_v2_potency",
        "entry": "predict.py",
        "weights": "kfm_v2_lsl_20260810.joblib",
        "ram_gb": 5.0,
        "blurb": "Which of two compounds is more potent against one kinase",
    },
    "selectivity": {
        "bundle": "SeqALigSeqB_v2_selectivity",
        "entry": "predict.py",
        "weights": "kfm_v2_sls_20260810.joblib",
        "ram_gb": 22.0,
        "blurb": "Which of two kinases binds one compound more tightly",
    },
}


MODELS_DIRNAME = "kfm-models"


def cache_root() -> str:
    """Where bundles live. See the module docstring.

    KFM_HOME wins outright. Otherwise walk up from the current directory looking
    for an existing `kfm-models`, so a command run from a subdirectory of the
    project still finds the models downloaded at its root. If none exists
    anywhere above, a fresh download goes to `./kfm-models`.
    """
    env = os.environ.get("KFM_HOME")
    if env:
        return os.path.abspath(os.path.expanduser(env))

    here = os.path.abspath(os.getcwd())
    node = here
    while True:
        candidate = os.path.join(node, MODELS_DIRNAME)
        if os.path.isdir(candidate):
            return candidate
        parent = os.path.dirname(node)
        if parent == node:            # reached the filesystem root
            break
        node = parent
    return os.path.join(here, MODELS_DIRNAME)


def bundle_dir(model: str) -> str:
    """Where this model's bundle lives, wherever it came from.

    Resolution order, first hit wins:
      1. KFM_BUNDLE_<MODEL>   — an explicit directory, for testing against a
                                bundle you already have
      2. <cache_root>/<bundle-name>
    """
    spec = MODELS[model]
    override = os.environ.get("KFM_BUNDLE_" + model.upper())
    if override:
        return os.path.abspath(os.path.expanduser(override))
    return os.path.join(cache_root(), spec["bundle"])


def is_present(model: str) -> bool:
    """Is this bundle usable: its own predict module, and a forest to load.

    The weights are matched by EXTENSION, not by an exact filename. Each model's
    file is named for what it is -- version, architecture, build date, e.g.
    kfm_v2_lsl_20260810.joblib -- so pinning the name here would mean this
    package had to be upgraded in lockstep with every model release, and an
    older copy would refuse a perfectly good bundle. `spec["weights"]` remains
    as the name a fresh download writes and as the legacy name still used by the
    deployed services.
    """
    import glob
    d = bundle_dir(model)
    spec = MODELS[model]
    if not os.path.isfile(os.path.join(d, spec["entry"])):
        return False
    return bool(glob.glob(os.path.join(d, "*.joblib")))


def total_ram_gb():
    """Physical RAM in GB, or None if it cannot be determined.

    Used only to warn. A wrong guess must never stop someone running the model:
    containers and VMs routinely misreport, and refusing to load on a machine
    that would in fact have coped is worse than letting it try.
    """
    try:
        if hasattr(os, "sysconf") and "SC_PAGE_SIZE" in os.sysconf_names \
                and "SC_PHYS_PAGES" in os.sysconf_names:
            return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3)
    except (ValueError, OSError):
        pass
    try:                                             # macOS
        import subprocess
        out = subprocess.run(["sysctl", "-n", "hw.memsize"],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return int(out.stdout.strip()) / (1024 ** 3)
    except Exception:                                # noqa: BLE001
        pass
    return None


def check_ram(model: str) -> None:
    """Warn -- loudly, before a long load -- if this machine looks too small.

    These forests expand roughly 7x from their file size when unpickled: the
    selectivity model is 2.8 GB on disk and 22 GB resident. A machine without
    the headroom does not fail cleanly, it thrashes swap for several minutes and
    is then killed by the OS, which reads as "the model is broken" rather than
    "this laptop is too small".
    """
    need = MODELS[model].get("ram_gb")
    have = total_ram_gb()
    if not need or have is None:
        return
    if have < need * 1.15:                # a little headroom for the interpreter
        sys.stderr.write(
            "\n  WARNING: the %s model needs about %.0f GB of RAM once loaded "
            "and this\n           machine has %.0f GB. It will likely swap "
            "heavily or be killed.\n           The file is only %.1f GB on "
            "disk; these forests expand about 7x\n           when they are "
            "loaded.\n\n" % (model, need, have, need / 7.0))


def _sha256(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _fetch(url: str, dest: str, expect_sha: str | None = None) -> None:
    """Download with a progress line, then verify.

    Downloads to a .part file and renames only on success, so an interrupted
    transfer can never leave a half-written model that loads and scores wrongly.
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r:
        total = int(r.headers.get("Content-Length") or 0)
        got = 0
        with open(tmp, "wb") as out:
            while True:
                buf = r.read(1 << 20)
                if not buf:
                    break
                out.write(buf)
                got += len(buf)
                if total:
                    pct = 100.0 * got / total
                    sys.stderr.write(
                        f"\r  {os.path.basename(dest)}  {got/1e6:,.0f} / "
                        f"{total/1e6:,.0f} MB  ({pct:4.1f}%)")
                else:
                    sys.stderr.write(
                        f"\r  {os.path.basename(dest)}  {got/1e6:,.0f} MB")
                sys.stderr.flush()
    sys.stderr.write("\n")
    if expect_sha:
        actual = _sha256(tmp)
        if actual != expect_sha:
            os.remove(tmp)
            raise RuntimeError(
                f"checksum mismatch for {os.path.basename(dest)}\n"
                f"  expected {expect_sha}\n  actual   {actual}\n"
                "The download is corrupt; nothing was installed.")
    os.replace(tmp, dest)


def download(model: str, accept_licence: bool = False) -> str:
    """Fetch a bundle into the cache. Returns its directory."""
    spec = MODELS[model]
    dest = bundle_dir(model)

    if not accept_licence:
        raise SystemExit(
            "\nThe model weights are licensed separately from this code.\n\n"
            "  RESEARCH AND EVALUATION USE ONLY.\n"
            "  Free for academic and non-profit research and teaching.\n\n"
            "  TWO LIMITS APPLY TO EVERYONE, INCLUDING ACADEMIC USERS:\n"
            "    * No licence is granted to the Eidogen-Sertanty Kinase\n"
            "      Knowledgebase itself, to any part of it, or to any data in\n"
            "      it. Running a model trained on it is not permission to use\n"
            "      it; that is a separate commercial agreement.\n"
            "    * No reverse engineering, and no extracting, reconstructing\n"
            "      or approximating the Knowledgebase from these weights or\n"
            "      their outputs, including by bulk querying.\n\n"
            "Full terms: LICENSE-MODELS.txt\n\n"
            "Re-run with --accept-licence to accept and download.\n")

    manifest_url = f"{BASE_URL}/{spec['bundle']}/manifest.json"
    print(f"Fetching {spec['bundle']} into {dest}")
    with urllib.request.urlopen(manifest_url, timeout=60) as r:
        manifest = json.load(r)

    os.makedirs(dest, exist_ok=True)
    for item in manifest["files"]:
        target = os.path.join(dest, item["name"])
        if os.path.isfile(target) and item.get("sha256"):
            if _sha256(target) == item["sha256"]:
                print(f"  {item['name']} already present, verified")
                continue
        _fetch(f"{BASE_URL}/{spec['bundle']}/{item['name']}", target,
               item.get("sha256"))
    print(f"\n{spec['bundle']} ready in {dest}")
    return dest


# The versions these forests were fitted with. A joblib-pickled random forest is
# a pickled object graph: loaded by a different scikit-learn it either dies with
# an unreadable traceback or, far worse, loads and scores differently.
PINNED = {"scikit-learn": "1.7.2", "numpy": "2.2.6", "joblib": "1.5.3"}


def check_libraries() -> None:
    """Fail clearly if this interpreter is not the one the install built.

    The failure this prevents is specific and common: `./install.sh` builds a
    correct environment in ./env, and then the user types `python -m kfm` in a
    shell whose `python` is a base conda install with scikit-learn 1.2 and
    numpy 1.x. What they got was `ModuleNotFoundError: No module named
    'numpy._core'` out of the depths of pickle -- which names neither the real
    problem nor the fix.

    Set KFM_SKIP_VERSION_CHECK=1 to bypass this deliberately.
    """
    if os.environ.get("KFM_SKIP_VERSION_CHECK"):
        return
    try:
        import importlib.metadata as md
        found = {name: md.version(name) for name in PINNED}
    except Exception:                                        # noqa: BLE001
        return                                               # never block on this
    wrong = {n: (v, PINNED[n]) for n, v in found.items() if v != PINNED[n]}
    if not wrong:
        return
    lines = [
        "",
        "This Python is not the one the install built, and its libraries do not",
        "match the versions these models were fitted with:",
        "",
    ]
    for n, (got, want) in sorted(wrong.items()):
        lines.append(f"    {n:<14} found {got:<10} needs {want}")
    lines += [
        "",
        f"    interpreter    {sys.executable}",
        "",
        "Loading a model with the wrong scikit-learn can score DIFFERENTLY rather",
        "than failing, so this stops here instead of returning a number.",
        "",
        "Run it through the wrapper, which always uses the right environment:",
        "",
        "    ./kfm.sh potency -t ABL1 -l @examples/ABL1_potency.smi",
        "",
        "or activate the environment first:",
        "",
        "    source ./env/bin/activate       # conda users: conda activate ./env",
        "",
        "If you have not built one yet:   ./install.sh potency",
        "",
    ]
    raise SystemExit("\n".join(lines))


def load_predict(model: str):
    """Import the BUNDLE'S OWN predict module and return it.

    The bundle directory goes on sys.path because the selectivity bundle's
    `embed.py` imports `predict` by name; loading the file in isolation would
    break that import.
    """
    if not is_present(model):
        raise SystemExit(
            f"The {model} model is not downloaded yet.\n"
            f"  Looked in: {bundle_dir(model)}\n\n"
            f"Get it with:  kfm download {model} --accept-licence\n"
            f"Or point at a bundle you already have:\n"
            f"  export KFM_BUNDLE_{model.upper()}=/path/to/bundle\n")
    check_libraries()
    check_ram(model)
    d = bundle_dir(model)
    if d not in sys.path:
        sys.path.insert(0, d)
    name = f"kfm_bundle_{model}"
    if name in sys.modules:
        return sys.modules[name]
    spec_ = importlib.util.spec_from_file_location(
        name, os.path.join(d, MODELS[model]["entry"]))
    mod = importlib.util.module_from_spec(spec_)
    sys.modules[name] = mod
    spec_.loader.exec_module(mod)
    return mod
