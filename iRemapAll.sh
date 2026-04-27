#!/bin/bash
#SBATCH --job-name=iRemapAll
#SBATCH --output=/PHShome/kv1026/p05/logs/iRemapAll.out
#SBATCH --error=/PHShome/kv1026/p05/logs/iRemapAll.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krisna_vandyke@mail.dfci.harvard.edu
#SBATCH --export=ALL

set -euo pipefail

echo "[$(date)] Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Shell: $SHELL"

# 0) Ensure Lmod (module) is initialized even in clean shells
if ! command -v module >/dev/null 2>&1; then
  # Try common init locations (safe to source if present)
  [ -f /etc/profile.d/lmod.sh ]    && . /etc/profile.d/lmod.sh
  [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
  [ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash
fi

# Sanity check
if ! command -v module >/dev/null 2>&1; then
  echo "FATAL: Lmod not initialized (no 'module' command). Check system profiles." >&2
  exit 2
fi

# 1) Clean up env that can confuse R
unset R_HOME || true

# 2) Try to load R (ignore cache so stale indexes don't bite us)
echo "Initializing MODULEPATH:"
echo "${MODULEPATH:-<empty>}" | tr ':' '\n' || true

module --ignore-cache load R/4.4.0 || {
  echo "WARN: 'module load R/4.4.0' failed. Showing spider info for R:" >&2
  module spider R || true
  echo "FATAL: Could not load an R module. Adjust the version to one listed above." >&2
  exit 3
}

# 3) Log what we actually got
echo "Loaded modules:"
module list 2>&1 || true

echo "Rscript location: $(command -v Rscript || echo 'not found')"
Rscript --version || { echo "FATAL: Rscript still not found after module load." >&2; exit 4; }

# 4) Run your script
echo "[$(date)] Starting R job…"
Rscript /PHShome/kv1026/p05/iRemapAll.R
echo "[$(date)] Done."
