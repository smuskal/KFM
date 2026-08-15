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
  echo "Could not find conda's profile hook. Either install miniforge/miniconda," >&2
  echo "or use the venv path instead, which needs no conda at all:" >&2
  echo "    KFM_FORCE_VENV=1 ./install.sh potency" >&2
  exit 1
fi
echo "conda   : $FOUND_HOOK"

mkdir -p "$SCRATCH"
cd "$SCRATCH"

# ---------------------------------------------------------------------------
echo
echo "=== 1/7  environment, created INSIDE the scratch directory ==="
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
echo "=== 2/7  weights, downloaded INSIDE the scratch directory ==="
export KFM_HOME="$SCRATCH/models"
# Any KFM_BUNDLE_* left in the caller's shell would silently point this test at
# a bundle somewhere else, which would make the whole exercise meaningless.
unset KFM_BUNDLE_POTENCY KFM_BUNDLE_SELECTIVITY 2>/dev/null || true
echo "  KFM_HOME=$KFM_HOME"
cd "$REPO"
"$PY" -m kfm download potency --accept-licence

# ---------------------------------------------------------------------------
echo
echo "=== 3/7  the worked example from the potency report ==="
echo "        bosutinib must come first, with score 0.847"
"$PY" -m kfm potency --target ABL1 \
  -l "COc1cc(Nc2c(cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)C#N)c(Cl)cc1Cl bosutinib" \
  -l "CC(C)n1nc(c2cccnc2)c3c(N)ncnc13 PP1-type"

# ---------------------------------------------------------------------------
echo
echo "=== 4/7  the test suite ==="
"$PY" -m pytest "$REPO/tests" -q

# ---------------------------------------------------------------------------
# The two add-ons, driven end to end. The suite covers them too, but only when a
# bundle is installed -- and this script is the one place we know one is. Both
# write into the scratch directory and neither touches ./kfm-models.
echo
echo "=== 5/7  kfm buildnew: fit on your data alone ==="
# Potency only. With no --recipe, buildnew reads the feature recipe from the
# bundle matching --layout, and this script deliberately downloads potency only
# -- selectivity is 2.84 GB and needs 22 GB of RAM to load. The suite covers the
# selectivity layout for anyone who has that bundle installed.
for layout in potency; do
  "$PY" -m kfm buildnew --layout "$layout" \
    --data "$REPO/examples/measurements_example.csv" \
    --out "$SCRATCH/buildnew_$layout" --trees 8 --allow-small
  "$PY" - "$SCRATCH/buildnew_$layout" "$layout" <<'PYCODE'
import json, os, sys
out, layout = sys.argv[1], sys.argv[2]
man = json.load(open(os.path.join(out, "MANIFEST.json")))
assert man["layout"] == layout, man.get("layout")
assert "derived_from" not in man, "buildnew must not claim KFM weights"
assert any(f.endswith(".joblib") for f in os.listdir(out)), "no model written"
print(f"    OK  {layout}: model written, no KFM weights claimed")
PYCODE
done

# ---------------------------------------------------------------------------
echo
echo "=== 6/7  kfm extend: merge your data into the released model ==="
"$PY" -m kfm extend --model potency \
  --data "$REPO/examples/extend_potency_example.csv" \
  --out "$SCRATCH/extended" --trees 5 --allow-small --label "contained-test"
"$PY" - "$SCRATCH/extended" <<'PYCODE'
import json, os, sys
out = sys.argv[1]
man = json.load(open(os.path.join(out, "MANIFEST.json")))
rec = man["extensions"][-1]
assert rec["trees_added"] == 5, rec["trees_added"]
assert rec["usable_comparisons"] == 16, rec["usable_comparisons"]
assert "EXTENDED_MODEL_WARNING" in man, "released figures were not withdrawn"
assert not os.path.exists(os.path.join(out, "reference_predictions.json"))
print(f"    OK  merged to {rec['n_estimators_after']} trees, "
      f"figures withdrawn, references retired")
PYCODE

# ---------------------------------------------------------------------------
echo
echo "=== 7/7  containment check ==="
echo "  everything created by this run:"
du -sh "$SCRATCH/env" "$SCRATCH/models" \
      "$SCRATCH"/buildnew_* "$SCRATCH/extended" 2>/dev/null | sed 's/^/    /'
echo
echo "  nothing was written to the default cache:"
if [ -e "$HOME/.cache/kfm" ]; then
  echo "    ~/.cache/kfm EXISTS - left by a version before models moved to ./kfm-models"
else
  echo "    ~/.cache/kfm does not exist (correct: models now live in ./kfm-models)"
fi
echo
echo "  to remove every trace of this test:"
echo "    rm -rf $SCRATCH"
