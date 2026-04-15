#!/usr/bin/env bash
set -euo pipefail

# Probe whether mature chain-building toolchains are available on CSC.
# This script is read-only (no file generation outside a temporary test).

echo "=== CSC liftover toolchain probe ==="
echo "Time: $(date)"
echo

module load minimap2 2>/dev/null || true
module load ucsc-tools 2>/dev/null || true
module load samtools 2>/dev/null || true
module load lastz 2>/dev/null || true
module load kentutils 2>/dev/null || true
module load bioinfo-tools 2>/dev/null || true

check_cmd() {
    local cmd="$1"
    if command -v "${cmd}" >/dev/null 2>&1; then
        echo "[OK] ${cmd}: $(command -v "${cmd}")"
    else
        echo "[--] ${cmd}: not found"
    fi
}

echo "--- Core commands ---"
check_cmd minimap2
check_cmd paftools.js
check_cmd liftOver
check_cmd axtChain
check_cmd chainSort
check_cmd chainNet
check_cmd netChainSubset
check_cmd lastz
check_cmd samtools
echo

echo "--- paftools.js chain support ---"
if command -v paftools.js >/dev/null 2>&1; then
    if paftools.js 2>&1 | grep -q -w "chain"; then
        echo "[OK] paftools.js supports 'chain' subcommand"
    else
        echo "[--] paftools.js present but no 'chain' subcommand"
    fi
else
    echo "[--] paftools.js not found"
fi
echo

echo "--- Suggested interpretation ---"
if command -v paftools.js >/dev/null 2>&1 && paftools.js 2>&1 | grep -q -w "chain"; then
    echo "[RECOMMENDED] You can use minimap2 + paftools.js chain directly."
elif command -v axtChain >/dev/null 2>&1 && command -v lastz >/dev/null 2>&1; then
    echo "[POSSIBLE] UCSC-style chain workflow may be feasible (lastz + axtChain + chainNet)."
else
    echo "[BLOCKED] No confirmed mature chain generation path detected."
    echo "          Try a different module version, conda env, or container with paftools chain."
fi

