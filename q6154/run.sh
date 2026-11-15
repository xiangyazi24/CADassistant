#!/usr/bin/env bash
set -euo pipefail
mkdir -p q6154
{
  date -u
  uname -a
  python3 --version
  python3 - <<'PY'
print('PUBLIC_RUNNER_OK')
PY
} > q6154/probe.txt
cat q6154/probe.txt
