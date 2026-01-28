#!/usr/bin/env bash
set -euo pipefail

# Probe available module versions on Anvil and test an Intel-based stack.
if ! command -v module >/dev/null 2>&1; then
  # Lmod typically provides module through shell initialization; try to source it.
  if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
fi

module purge

echo "==> Available Intel modules"
module spider intel

echo "==> Loading Intel (default)"
module load intel
module list

echo "==> Available MPI modules (after Intel)"
module spider openmpi
module spider impi

echo "==> Available IO libraries (after Intel)"
module spider hdf5
module spider netcdf
module spider pnetcdf

echo "==> Available build tools"
module spider cmake
