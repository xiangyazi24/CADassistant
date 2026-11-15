#!/usr/bin/env bash
set -euo pipefail
mkdir -p q6154
python3 -m pip install --quiet sympy
cp q6154/phase1.py q6154/phase1-fast.py
sed -i 's/^N = 5000$/N = 700/; s/^FIT = 240$/FIT = 180/' q6154/phase1-fast.py
python3 q6154/phase1-fast.py
