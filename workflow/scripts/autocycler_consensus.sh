#!/usr/bin/env bash
set -euo pipefail
work=$1
candidates=$2
reads=$3
threads=$4
output=$5
shopt -s nullglob
assemblies=("$candidates"/assemblies/*.fasta)
if (( ${#assemblies[@]} < 2 )); then
    echo "Autocycler requires at least two successful input assemblies." >&2
    exit 1
fi
mkdir -p "$work"
attempt=$(mktemp -d "$work/attempt_XXXXXX")
echo "Autocycler working directory: $attempt"
graph="$attempt/autocycler_out"
autocycler compress -i "$candidates/assemblies" -a "$graph" -t "$threads"
autocycler cluster -a "$graph"
clusters=("$graph"/clustering/qc_pass/cluster_*)
if (( ${#clusters[@]} == 0 )); then
    echo "No clusters passed Autocycler QC." >&2
    exit 1
fi
for cluster in "${clusters[@]}"; do
    autocycler trim -c "$cluster" -t "$threads"
    autocycler resolve -c "$cluster"
done
autocycler combine -a "$graph" -i "$graph"/clustering/qc_pass/cluster_*/5_final.gfa -r "$reads" -t "$threads"
test -s "$graph/consensus_assembly.fasta"
cp "$graph/consensus_assembly.fasta" "$output"
