#!/usr/bin/env bash
#
#   ./kfm.sh potency -t ABL1 -l @examples/ABL1_potency.smi
#
# Runs the tool using the environment ./install.sh built, whatever shell you are
# in and whether or not you have activated anything.
#
# This exists because the most likely way to get a wrong answer is to type
# `python -m kfm` in a terminal whose `python` is some other installation. That
# python has its own scikit-learn, and a joblib-pickled forest loaded by the
# wrong scikit-learn either dies with an unreadable traceback or -- worse --
# loads and scores differently. Using this wrapper removes the choice.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -x "$HERE/env/bin/python" ]; then
  exec "$HERE/env/bin/python" -m kfm "$@"
fi

if [ -x "$HERE/.venv/bin/python" ]; then
  exec "$HERE/.venv/bin/python" -m kfm "$@"
fi

cat >&2 <<EOF
No environment found in this directory.

Build one first, which also downloads the models:

    cd "$HERE"
    ./install.sh potency

EOF
exit 1
