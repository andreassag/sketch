#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/sketch"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== 1. Sketching genomes ==="
$BIN sketch -i data/strain_A.fa -o strain_A.skh -k 21 -s 11
$BIN sketch -i data/strain_B.fa -o strain_B.skh -k 21 -s 11
$BIN sketch -i data/strain_C.fa -o strain_C.skh -k 21 -s 11

echo -e "\n=== 2. Pairwise Distance & ANI ==="
echo "--- Strain A vs Strain B (expect ~99% ANI) ---"
$BIN dist strain_A.skh strain_B.skh

echo -e "\n--- Strain A vs Strain C (expect ~90% ANI) ---"
$BIN dist strain_A.skh strain_C.skh

echo -e "\n=== 3. All-vs-All ANI Matrix ==="
$BIN matrix strain_A.skh strain_B.skh strain_C.skh

echo -e "\n=== 4. Neighbor-Joining Phylogenetic Tree ==="
$BIN tree strain_A.skh strain_B.skh strain_C.skh

rm -f strain_A.skh strain_B.skh strain_C.skh
echo -e "\nExample run completed successfully!"
