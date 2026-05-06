#!/bin/bash

# Binding script for Perlmutter GPU

set -euo pipefail

# Local rank on the node
lrank=$(( $SLURM_LOCALID % 4 ))

# GPU and NIC binding
export MPICH_OFI_NIC_POLICY="USER"
export MPICH_OFI_NIC_MAPPING="0:3;1:2;2:1;3:0"

# Print what each rank will run
echo "rank=$SLURM_PROCID localid=$SLURM_LOCALID lrank=$lrank cmd: $*" >&2

# Binding
case "$lrank" in
 0) exec numactl --physcpubind=0-15,64-79    --membind=0 "$@" ;;
 1) exec numactl --physcpubind=16-31,80-95   --membind=1 "$@" ;;
 2) exec numactl --physcpubind=32-47,96-111  --membind=2 "$@" ;;
 3) exec numactl --physcpubind=48-63,112-127 --membind=3 "$@" ;;
esac
