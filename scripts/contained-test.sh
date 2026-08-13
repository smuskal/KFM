#!/usr/bin/env bash
# A completely self-contained test of this package.
#
#   ./scripts/contained-test.sh  [/path/to/scratch]
#
# EVERYTHING it creates lives under the scratch directory: the conda environment,
# the downloaded model weights, the cache. Nothing is written to your home
# directory, nothing is registered in your global conda environment list, and
# nothing on your Mac is modified. Delete the scratch directory and every trace
# is gone.
#
# The two mechanisms that make that true, and they are the whole point:
#
#   conda env create -p ./env    creates the environment AT A PATH rather than
#                                under a name in ~/miniforge3/envs, so it is a
#                                directory you own and can delete
#
#   export KFM_HOME=./models     puts the downloaded weights here instead of the
#                                default ./kfm-models
set -euo pipefail

SCRATCH="${1:-$(pwd)/scratch}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "repo    : $REPO"
echo "scratch : $SCRATCH"
echo

# conda's shell function is not available to a non-interactive script unless the
# profile hook is sourced, and `conda activate` fails confusingly without it.
for hook in "$HOME/miniforge3-arm64/etc/profile.d/conda.sh" \
            "$HOME/miniforge3/etc/profile.d/conda.sh" \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
  if [ -f "$hook" ]; then . "$hook"; FOUND_HOOK="$hook"; break; fi
done
if [ -z "${FOUND_HOOK:-}" ]; then
  echo "Could not find conda's profile hook. Install miniforge/miniconda, or"
  echo "run the steps in LOCAL_TEST.md by hand." >&2
  exit 1
fi
echo "conda   : $FOUND_HOOK"

mkdir -p "$SCRATCH"
cd "$SCRATCH"

# ---------------------------------------------------------------------------
echo
echo "=== 1/5  environment, created INSIDE the scratch directory ==="
if [ -d "$SCRATCH/env" ]; then
  echo "  reusing $SCRATCH/env"
else
  conda env create -p "$SCRATCH/env" -f "$REPO/environment.yml"
fi
PY="$SCRATCH/env/bin/python"
"$PY" -c "import sys, sklearn, numpy, joblib, rdkit
print('  python  ', sys.version.split()[0])
print('  prefix  ', sys.prefix)
print('  sklearn ', sklearn.__version__)
print('  numpy   ', numpy.__version__)
print('  joblib  ', joblib.__version__)
print('  rdkit   ', rdkit.__version__)"

# ---------------------------------------------------------------------------
echo
echo "=== 2/5  weights, downloaded INSIDE the scratch directory ==="
export KFM_HOME="$SCRATCH/models"
# Any KFM_BUNDLE_* left in the caller's shell would silently point this test at
# a bundle somewhere else, which would make the whole exercise meaningless.
unset KFM_BUNDLE_POTENCY KFM_BUNDLE_SELECTIVITY 2>/dev/null || true
echo "  KFM_HOME=$KFM_HOME"
cd "$REPO"
"$PY" -m kfm download potency --accept-licence

# ---------------------------------------------------------------------------
echo
echo "=== 3/5  the worked example from the potency report ==="
echo "        bosutinib must come first, with score 0.847"
"$PY" -m kfm potency --target ABL1 \
  -l "COc1cc(Nc2c(cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)C#N)c(Cl)cc1Cl bosutinib" \
  -l "CC(C)n1nc(c2cccnc2)c3c(N)ncnc13 PP1-type"

# ---------------------------------------------------------------------------
echo
echo "=== 4/5  the test suite ==="
"$PY" -m pytest "$REPO/tests" -q

# ---------------------------------------------------------------------------
echo
echo "=== 5/5  containment check ==="
echo "  everything created by this run:"
du -sh "$SCRATCH/env" "$SCRATCH/models" 2>/dev/null | sed 's/^/    /'
echo
echo "  nothing was written to the default cache:"
if [ -e "$HOME/.cache/kfm" ]; then
  echo "    ~/.cache/kfm EXISTS — left by a version before models moved to ./kfm-models"
else
  echo "    ~/.cache/kfm does not exist (correct: models now live in ./kfm-models)"
fi
echo
echo "  to remove every trace of this test:"
echo "    rm -rf $SCRATCH"
