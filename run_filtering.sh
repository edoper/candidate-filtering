#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────
#  End-to-end clinical filtering with Pangolin splice scoring.
#
#  Pass 1: filtering_r.pl emits <proband>.pangolin_input.csv (rare, panel,
#          MANE, whitelisted "structural-pass" variants) for any proband that
#          does not yet have a Pangolin score map.
#  Pangolin: compute de-novo splice scores on those variants (GPU).
#  Parse:    build <proband>.pangolin.tsv (chr-pos-ref-alt -> max|delta|).
#  Pass 2:   filtering_r.pl produces <proband>.<panel>.candidatos using the
#            scores (pangolin_score column + splice rescue at >= 0.5).
#  Cleanup:  all Pangolin scratch is deleted; only <proband>.<panel>.candidatos
#            and the annotated VCFs (*.germline.vep.vcf.gz + .tbi + _summary.html)
#            survive. Pangolin is therefore recomputed every run (cheap — only
#            the few hundred structural-pass candidates).
# ──────────────────────────────────────────────────────────────────────────
set -euo pipefail
# Optional arg $1 = genes-of-interest file (one symbol per line) forwarded to
# filtering_r.pl; if omitted the default panel (g4e-2025.txt) is used.
# Working directory holds the *.germline.vep.vcf.gz data + reference files
# (defaults to this script's own directory; override with $WORKDIR).
GENES="${1:-}"
cd "${WORKDIR:-$(dirname "$(readlink -f "$0")")}"

# Optional: force specific sample(s) as proband, overriding filename-based
# auto-discovery.  e.g.  PROBAND="EPIC280M" bash run_filtering.sh
#                        PROBAND="EPIC280 EPIC280M" bash run_filtering.sh genes.txt
PROBAND_ARGS=()
for _p in ${PROBAND:-}; do PROBAND_ARGS+=(--proband "$_p"); done
FWD=("${PROBAND_ARGS[@]+"${PROBAND_ARGS[@]}"}" ${GENES:+"$GENES"})

# Configurable environment (defaults match the original setup).
source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
conda activate "${PANGOLIN_ENV:-pangolin}"

FA="${PANGOLIN_FASTA:-$HOME/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa}"
DB="${PANGOLIN_DB:-$HOME/vep_refs/pangolin/gencode.v38.annotation.db}"

echo "===== Pass 1: emit Pangolin candidate inputs ====="
perl filtering_r.pl "${FWD[@]}"

echo "===== Pangolin scoring ====="
shopt -s nullglob
for csv in *.pangolin_input.csv; do
    proband="${csv%.pangolin_input.csv}"
    tsv="$proband.pangolin.tsv"
    echo "[pangolin] scoring $proband ($(($(wc -l < "$csv") - 1)) variants) ..."
    pangolin "$csv" "$FA" "$DB" "$proband.pangolin" -c CHROM,POS,REF,ALT
    perl parse_pangolin.pl "$proband.pangolin.csv" > "$tsv"
    echo "[pangolin] -> $tsv ($(wc -l < "$tsv") variants scored)"
done

echo "===== Pass 2: final filtering with splice scores ====="
perl filtering_r.pl "${FWD[@]}"

# Keep only the final tables + annotated VCFs; drop all regenerable Pangolin
# scratch. Runs only on success (set -e aborts earlier), so a failed run leaves
# intermediates in place for debugging.
echo "===== Cleanup: removing Pangolin intermediates ====="
rm -f -- *.pangolin_input.csv *.pangolin.csv *.pangolin.tsv *.pangolin.md5

echo "===== done ====="
