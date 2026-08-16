#!/usr/bin/env bash
#
#   ./install.sh                 # everything: environment + both version 2 models
#   ./install.sh potency         # environment + just the potency model
#   ./install.sh all             # environment + version 2 AND the legacy version 1
#
# One command does the whole job. EVERYTHING lands inside this directory: the
# environment in ./env, the model weights in ./kfm-models. Nothing is written to
# your home directory and no environment you already have is modified. Delete
# this directory and every trace is gone.
#
# It prefers conda and falls back to python's own venv, so a missing or broken
# conda is not a dead end. It refuses to finish unless each model it installed
# returns its published answer - "it installed" and "it works" are the same
# statement here, not two different ones.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
ENVDIR="$HERE/env"

say()  { printf '\n\033[1m%s\033[0m\n' "$*"; }
ok()   { printf '\033[32m%s\033[0m\n' "$*"; }
fail() { printf '\n\033[31mFAILED: %s\033[0m\n' "$*" >&2; exit 1; }

# --- which models -----------------------------------------------------------
# Default to the current release, both models, because installing one and
# discovering later that the other is missing is the most common way to waste
# an afternoon. Names are validated BEFORE anything is downloaded: a typo that
# surfaces after a 3 GB transfer is a bad way to learn you meant "selectivity".
if [ $# -eq 0 ]; then
  MODELS=(potency selectivity)
elif [ "${1:-}" = "all" ]; then
  MODELS=(potency selectivity v1)
else
  MODELS=("$@")
fi
for m in "${MODELS[@]}"; do
  case "$m" in
    potency|selectivity|v1) ;;
    *) fail "unknown model '$m'.
        Valid names are: potency, selectivity, v1
        Or use:  ./install.sh          (both version 2 models)
                 ./install.sh all      (version 2 and version 1)" ;;
  esac
done

say "Kinase Foundation Model - self-contained install"
echo "  installing:    ${MODELS[*]}"
echo "  everything in: $HERE"
echo "  environment:   $ENVDIR"
echo "  model weights: $HERE/kfm-models"

# --- the licence, once, for everything --------------------------------------
cat <<'EOF'

  The model weights are licensed separately from this code.

    RESEARCH AND EVALUATION USE ONLY. Free for academic and non-profit
    research and teaching.

    No licence is granted to the Eidogen-Sertanty Kinase Knowledgebase
    itself, and no reverse engineering or bulk extraction of it from these
    weights or their outputs. Full terms: LICENSE-MODELS.txt

  Continuing this install accepts those terms for every model above.
EOF

# --- how much room is needed ------------------------------------------------
# Checked before downloading, not after. Selectivity is 2.84 GB to fetch and
# 22 GB to load, and a machine that cannot run it should learn that now.
NEED_RAM=0
for m in "${MODELS[@]}"; do
  case "$m" in
    selectivity) [ "$NEED_RAM" -lt 22 ] && NEED_RAM=22 ;;
    v1)          [ "$NEED_RAM" -lt 7 ]  && NEED_RAM=7  ;;
    potency)     [ "$NEED_RAM" -lt 5 ]  && NEED_RAM=5  ;;
  esac
done
HAVE_RAM=$(
  if [ "$(uname)" = "Darwin" ]; then sysctl -n hw.memsize 2>/dev/null | awk '{printf "%d", $1/1073741824}'
  else awk '/MemTotal/ {printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null; fi
) || HAVE_RAM=0
if [ "${HAVE_RAM:-0}" -gt 0 ] && [ "$HAVE_RAM" -lt "$NEED_RAM" ]; then
  printf '\n\033[33m  WARNING: %s needs about %s GB of RAM to load and this machine has %s GB.\n' \
         "${MODELS[*]}" "$NEED_RAM" "$HAVE_RAM"
  printf '  The download will work; loading the model will likely be killed by the OS.\n'
  printf '  Consider installing only the potency model:  ./install.sh potency\033[0m\n'
fi

# --- 1. the environment ------------------------------------------------------
# Created AT A PATH (-p), never by name: a named environment goes into your
# conda installation and stays there. A path environment is a directory you own.
if [ -x "$ENVDIR/bin/python" ]; then
  say "1/4  environment already exists, reusing it"
else
  CONDA=""
  if [ -z "${KFM_FORCE_VENV:-}" ]; then
    for hook in "$HOME/miniforge3/etc/profile.d/conda.sh" \
                "$HOME/miniforge3-arm64/etc/profile.d/conda.sh" \
                "$HOME/miniconda3/etc/profile.d/conda.sh" \
                "$HOME/anaconda3/etc/profile.d/conda.sh" \
                "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"; do
      [ -f "$hook" ] && { CONDA="$hook"; break; }
    done
  fi

  if [ -n "$CONDA" ]; then
    say "1/4  creating a conda environment in ./env"
    # shellcheck disable=SC1090
    . "$CONDA"
    conda env create -p "$ENVDIR" -f "$HERE/environment.yml" \
      || fail "conda could not build the environment. Delete ./env and re-run,
        or re-run as:  KFM_FORCE_VENV=1 ./install.sh ${MODELS[*]}"
  else
    say "1/4  using python's built-in venv in ./env"
    PY=""
    for c in python3.12 python3.11 python3.10 python3; do
      command -v "$c" >/dev/null 2>&1 && { PY="$c"; break; }
    done
    [ -n "$PY" ] || fail "no python 3.10-3.12 found on PATH.
        Install Python from https://www.python.org/downloads/ and open a new terminal."
    "$PY" -m venv "$ENVDIR" || fail "could not create a virtual environment."
    "$ENVDIR/bin/pip" install --quiet --upgrade pip
    "$ENVDIR/bin/pip" install --quiet -r "$HERE/requirements.txt" \
      || fail "pip could not install the pinned libraries."
  fi
fi

PY="$ENVDIR/bin/python"
[ -x "$PY" ] || fail "no interpreter at $PY"

# Not fatal, and deliberately so: `./kfm.sh` and `python -m kfm` both work from
# this directory whether or not the editable install lands. Only the bare `kfm`
# console script needs it. But swallowing the output entirely meant a failure
# here was completely invisible, so say so and carry on rather than hiding it.
if ! "$PY" -m pip install --quiet -e "$HERE" >/tmp/kfm_editable_install.log 2>&1; then
  printf '\n\033[33m  NOTE: the editable install did not complete.\n'
  printf '  ./kfm.sh and ./env/bin/python -m kfm still work; only the bare `kfm`\n'
  printf '  command is unavailable. Details: /tmp/kfm_editable_install.log\033[0m\n'
fi

# --- 2. what actually got installed -----------------------------------------
say "2/4  libraries installed"
"$PY" - <<'PYCODE'
import platform, sys, sklearn, numpy, joblib, rdkit
print(f"  python       {sys.version.split()[0]}  ({platform.machine()})")
print(f"  scikit-learn {sklearn.__version__}")
print(f"  numpy        {numpy.__version__}")
print(f"  joblib       {joblib.__version__}")
print(f"  rdkit        {rdkit.__version__}")
PYCODE

# --- 3. the models, one after another ---------------------------------------
say "3/4  downloading ${#MODELS[@]} model(s) into ./kfm-models"
for m in "${MODELS[@]}"; do
  echo
  echo "  --- $m ---"
  ( cd "$HERE" && "$PY" -m kfm download "$m" --accept-licence ) \
    || fail "the $m download did not complete. Re-run the same command: files that
        finished are checksummed and skipped, so it resumes rather than restarting."
done

# --- 4. prove each one works -------------------------------------------------
# An install that reports success without scoring anything is how a broken
# environment reaches a user. Each model must reproduce its published answer.
say "4/4  checking each model against its published answer"
BOS="COc1cc(Nc2c(cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)C#N)c(Cl)cc1Cl"
PP1="CC(C)n1nc(c2cccnc2)c3c(N)ncnc13"
DAS="Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1"

for m in "${MODELS[@]}"; do
  case "$m" in
    potency)
      OUT=$( cd "$HERE" && "$PY" -m kfm potency -t ABL1 -l "$BOS bosutinib" -l "$PP1 PP1-type" 2>&1 )
      echo "$OUT" | grep -q "0.847" \
        || { echo "$OUT" | head -6; fail "potency: on ABL1, bosutinib did not rank above the PP1-type compound at 0.847."; }
      ok "  potency      on ABL1, bosutinib over PP1-type 0.847 - matches the published report" ;;
    selectivity)
      OUT=$( cd "$HERE" && "$PY" -m kfm selectivity -t ABL1 -t GSK3B -l "$DAS dasatinib" 2>&1 )
      echo "$OUT" | grep -q "0.97" \
        || { echo "$OUT" | head -6; fail "selectivity: dasatinib did not prefer ABL1 at 0.97 over GSK3B."; }
      ok "  selectivity  dasatinib prefers ABL1 0.97 over GSK3B - matches the report" ;;
    v1)
      OUT=$( cd "$HERE" && "$PY" -m kfm v1 -t ABL1 -l "$BOS bosutinib" -l "$PP1 PP1-type" 2>&1 )
      echo "$OUT" | grep -q "8.69" \
        || { echo "$OUT" | head -6; fail "v1: bosutinib did not score 8.69 against ABL1."; }
      ok "  v1           bosutinib 8.69 vs ABL1 - matches the published report" ;;
  esac
done

say "Done. To use it:"
cat <<EOF
  cd $HERE
  ./kfm.sh potency -t ABL1 -l @examples/ABL1_potency.smi

Always use ./kfm.sh (or ./env/bin/python -m kfm). It runs the environment this
script just built. Typing plain \`python -m kfm\` or \`kfm\` uses whatever python
your shell happens to have - usually a base conda install with a different
scikit-learn, which cannot load these models.

  ./kfm.sh where                     # where the models are, and the RAM each needs
  ./kfm.sh selectivity -t ABL1 -t GSK3B -l "CCO ethanol"

To check the install yourself:
  ./env/bin/python -m pytest tests/ -q

To build on your own data (both run locally and transmit nothing):
  ./kfm.sh buildnew --help           # fit on your data alone
  ./kfm.sh extend   --help           # merge your data into a released model

Both take gene symbols with what is already installed. Only raw amino-acid
sequences for kinases the model has never seen need anything more:
  ./env/bin/pip install -r requirements-sequences.txt

To remove everything, delete this directory. Nothing else was touched.
EOF
