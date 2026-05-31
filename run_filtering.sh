#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────
#  End-to-end clinical filtering with Pangolin splice scoring.
#
#  Pass 1: filtering_r.pl emits <proband>.pangolin_input.csv (rare, panel,
#          MANE, whitelisted "structural-pass" variants) for any proband that
#          does not yet have a Pangolin score map.
#  Pangolin: compute de-novo splice scores on those variants (GPU).
#  Parse:    build <proband>.pangolin.tsv (chr-pos-ref-alt -> max|delta|).
#  Pass 2:   filtering_r.pl produces <proband>.germline.candidatos using the
#            scores (pangolin_score column + splice rescue at >= 0.5).
# ──────────────────────────────────────────────────────────────────────────
set -euo pipefail
# Run from the directory that holds the *.germline.vep.vcf.gz data + reference
# files (defaults to this script's own directory; override with arg $1).
cd "${1:-$(dirname "$(readlink -f "$0")")}"

# Configurable environment (defaults match the original setup).
source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
conda activate "${PANGOLIN_ENV:-pangolin}"

FA="${PANGOLIN_FASTA:-$HOME/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa}"
DB="${PANGOLIN_DB:-$HOME/vep_refs/pangolin/gencode.v38.annotation.db}"

echo "===== Pass 1: emit Pangolin candidate inputs ====="
perl filtering_r.pl

echo "===== Pangolin scoring ====="
shopt -s nullglob
for csv in *.pangolin_input.csv; do
    proband="${csv%.pangolin_input.csv}"
    tsv="$proband.pangolin.tsv"
    [ -e "$tsv" ] && { echo "[skip] $tsv exists"; continue; }
    echo "[pangolin] scoring $proband ($(($(wc -l < "$csv") - 1)) variants) ..."
    pangolin "$csv" "$FA" "$DB" "$proband.pangolin" -c CHROM,POS,REF,ALT
    perl parse_pangolin.pl "$proband.pangolin.csv" > "$tsv"
    echo "[pangolin] -> $tsv ($(wc -l < "$tsv") variants scored)"
done

echo "===== Pass 2: final filtering with splice scores ====="
perl filtering_r.pl

echo "===== done ====="
