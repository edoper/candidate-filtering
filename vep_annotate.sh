#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────
#  VEP annotation pipeline with:
#    - All plugins: LoF (LOFTEE), REVEL, AlphaMissense, EVE, CADD (SNV + indel)
#    - Custom annotation: gnomAD v4.1 joint (PASS + EXOMES_FILTERED + GENOMES_FILTERED)
#                         exposing AC_joint, AN_joint, AF_joint, nhomalt_joint, FILTER
#    - Custom annotation: ClinVar (NCBI weekly VCF, chr-prefixed) exposing
#                         CLNSIG, CLNREVSTAT, CLNDN  (independent of cache CLIN_SIG)
#    - Multiallelic sites are split (bcftools norm -m-any) before VEP so that
#      downstream chr-pos-ref-alt keys are unambiguous.
#
#  Usage:  vep_annotate.sh <input.vcf[.gz]> <output.vcf.gz>
#
#  Notes:
#    - Auto-detects chr vs no-chr convention in input and normalizes to 'chr'
#      (required because the gnomAD_min custom VCF uses 'chr' prefix).
#    - Uses system Perl 5.34 (/usr/bin/perl) — DOES NOT activate vep-loftee
#      conda env (its Perl 5.32 is missing List::MoreUtils).
#    - LD_PRELOAD ~/htslib/libhts.so resolves runtime symbols for LOFTEE's
#      Bio::DB::BigFile.so which does not link to libhts.
# ──────────────────────────────────────────────────────────────────────────

set -euo pipefail

INPUT="${1:?Usage: $0 <input.vcf[.gz]> <output.vcf.gz>}"
OUTPUT="${2:?Usage: $0 <input.vcf[.gz]> <output.vcf.gz>}"

# ── Paths ──
VEP="$HOME/ensembl-vep/vep"
VEP_DATA="$HOME/vep_data"
VEP_REFS="$HOME/vep_refs"
GNOMAD_VCF="$VEP_REFS/gnomAD_min/gnomAD.joint.v4.1.mane.all.vcf.gz"
CLINVAR_VCF="$VEP_REFS/clinvar/clinvar.chr.vcf.gz"
CADD_SNV="$VEP_REFS/CADD/whole_genome_SNVs.tsv.gz"            # CADD GRCh38 v1.7
CADD_INDEL="$VEP_REFS/CADD/gnomad.genomes.r4.0.indel.tsv.gz" # CADD GRCh38 v1.7 indels

# ── Environment ──
export PERL5LIB="$HOME/.vep/Plugins/loftee:$HOME/perl5/lib/perl5:${PERL5LIB:-}"
export LD_PRELOAD="$HOME/htslib/libhts.so${LD_PRELOAD:+:$LD_PRELOAD}"

# ── Sanity checks ──
[[ -x "$VEP"        ]] || { echo "ERROR: VEP not found at $VEP" >&2; exit 1; }
[[ -s "$GNOMAD_VCF" ]] || { echo "ERROR: gnomAD custom VCF not found: $GNOMAD_VCF" >&2; exit 1; }
[[ -s "$CADD_SNV"   ]] || { echo "ERROR: CADD SNV file not found: $CADD_SNV" >&2; exit 1; }
[[ -s "$CADD_INDEL" ]] || { echo "ERROR: CADD indel file not found: $CADD_INDEL" >&2; exit 1; }
[[ -s "$CLINVAR_VCF" ]] || { echo "ERROR: ClinVar custom VCF not found: $CLINVAR_VCF" >&2; exit 1; }
[[ -s "$INPUT"      ]] || { echo "ERROR: input not found or empty: $INPUT" >&2; exit 1; }

if ! /usr/bin/perl -MList::MoreUtils -MBio::DB::HTS -MBio::DB::BigFile -e1 2>/dev/null; then
  echo "ERROR: system Perl missing required modules (List::MoreUtils / Bio::DB::HTS / Bio::DB::BigFile)" >&2
  exit 1
fi

# ── Detect chromosome naming convention in input ──
FIRST_CHR=$(set +o pipefail; bcftools query -f '%CHROM\n' "$INPUT" 2>/dev/null | awk 'NR==1{print; exit}')
[[ -z "$FIRST_CHR" ]] && { echo "ERROR: no variants in $INPUT" >&2; exit 1; }

CLEANUP=()
if [[ "$FIRST_CHR" == chr* ]]; then
  echo "[vep] Input uses 'chr' prefix — no normalization needed"
  VEP_INPUT="$INPUT"
else
  echo "[vep] Input lacks 'chr' prefix — normalizing to temp file"
  VEP_INPUT="$(mktemp -u --suffix=.vcf.gz)"
  CLEANUP+=("$VEP_INPUT")

  # GRCh38 contig lengths (gnomAD canonical) — added to header so VEP/bcftools
  # don't emit warnings about undefined contigs.
  read -r -d '' CONTIGS <<'EOF' || true
##contig=<ID=chr1,length=248956422,assembly=GRCh38>
##contig=<ID=chr2,length=242193529,assembly=GRCh38>
##contig=<ID=chr3,length=198295559,assembly=GRCh38>
##contig=<ID=chr4,length=190214555,assembly=GRCh38>
##contig=<ID=chr5,length=181538259,assembly=GRCh38>
##contig=<ID=chr6,length=170805979,assembly=GRCh38>
##contig=<ID=chr7,length=159345973,assembly=GRCh38>
##contig=<ID=chr8,length=145138636,assembly=GRCh38>
##contig=<ID=chr9,length=138394717,assembly=GRCh38>
##contig=<ID=chr10,length=133797422,assembly=GRCh38>
##contig=<ID=chr11,length=135086622,assembly=GRCh38>
##contig=<ID=chr12,length=133275309,assembly=GRCh38>
##contig=<ID=chr13,length=114364328,assembly=GRCh38>
##contig=<ID=chr14,length=107043718,assembly=GRCh38>
##contig=<ID=chr15,length=101991189,assembly=GRCh38>
##contig=<ID=chr16,length=90338345,assembly=GRCh38>
##contig=<ID=chr17,length=83257441,assembly=GRCh38>
##contig=<ID=chr18,length=80373285,assembly=GRCh38>
##contig=<ID=chr19,length=58617616,assembly=GRCh38>
##contig=<ID=chr20,length=64444167,assembly=GRCh38>
##contig=<ID=chr21,length=46709983,assembly=GRCh38>
##contig=<ID=chr22,length=50818468,assembly=GRCh38>
##contig=<ID=chrX,length=156040895,assembly=GRCh38>
##contig=<ID=chrY,length=57227415,assembly=GRCh38>
##contig=<ID=chrM,length=16569,assembly=GRCh38>
EOF

  # Inject contig headers + prepend 'chr' to data + drop any pre-existing contig lines
  bcftools view "$INPUT" \
    | awk -v contigs="$CONTIGS" '
        BEGIN { OFS="\t" }
        /^##contig=/                  { next }
        /^#CHROM/                     { print contigs; print; next }
        /^#/                          { print; next }
                                      { $1 = "chr"$1; print }' \
    | bgzip > "$VEP_INPUT"
fi

# ── Split multiallelic sites (one ALT per record) ──
# Keeps downstream chr-pos-ref-alt keys unambiguous (filtering_r.pl, inheritance,
# Pangolin input). Split only (-m-any), no --check-ref, to avoid dropping records.
echo "[vep] Splitting multiallelic sites (bcftools norm -m-any)"
NORM_INPUT="$(mktemp -u --suffix=.vcf.gz)"
CLEANUP+=("$NORM_INPUT")
bcftools norm -m-any "$VEP_INPUT" -Oz -o "$NORM_INPUT"
VEP_INPUT="$NORM_INPUT"

# ── Run VEP ──
"$VEP" \
  --input_file       "$VEP_INPUT" \
  --output_file      "$OUTPUT" \
  --vcf \
  --compress_output  bgzip \
  --force_overwrite \
  --cache \
  --offline \
  --dir_cache        "$VEP_DATA" \
  --dir_plugins      "$HOME/.vep/Plugins" \
  --assembly         GRCh38 \
  --species          homo_sapiens \
  --fork             12 \
  --buffer_size      5000 \
  --mane_select \
  --hgvs \
  --symbol \
  --canonical \
  --biotype \
  --numbers \
  --domains \
  --regulatory \
  --variant_class \
  --check_existing \
  --plugin LoF,loftee_path:"$HOME/.vep/Plugins/loftee",human_ancestor_fa:"$VEP_REFS/loftee/GRCh38/human_ancestor.fa.gz",conservation_file:"$VEP_REFS/loftee/GRCh38/loftee.sql",gerp_bigwig:"$VEP_REFS/loftee/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
  --plugin REVEL,"$VEP_REFS/REVEL/new_tabbed_revel_grch38.tsv.gz" \
  --plugin AlphaMissense,file="$VEP_REFS/AlphaMissense/AlphaMissense_hg38.tsv.gz" \
  --plugin EVE,file="$VEP_REFS/EVE/eve_merged.vcf.gz" \
  --plugin CADD,"$CADD_SNV","$CADD_INDEL" \
  --custom file="$GNOMAD_VCF",short_name=gnomADmin,format=vcf,type=exact,fields=AC_joint%AN_joint%AF_joint%nhomalt_joint%FILTER \
  --custom file="$CLINVAR_VCF",short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNREVSTAT%CLNDN

# ── Index output ──
bcftools index -ft "$OUTPUT" 2>/dev/null || true

# ── Cleanup temp ──
for f in "${CLEANUP[@]}"; do rm -f "$f"; done

# ── Summary ──
# gnomADmin fields are *inside* CSQ (pipe-separated), not separate INFO tags.
# A hit is detected by the locator field "chr<n>:<pos>-<pos>" injected by VEP --custom.
N=$(bcftools view "$OUTPUT" 2>/dev/null | grep -vc "^#" || true)
GNOMAD_HITS=$(bcftools view "$OUTPUT" 2>/dev/null \
  | awk '/^#/{next} {print}' \
  | grep -cE 'chr[0-9XYM]+:[0-9]+-[0-9]+\|[0-9]+\|[0-9]+\|' || true)
SIZE=$(du -h "$OUTPUT" | cut -f1)

echo ""
echo "════════════════════════════════════════════════════"
echo "[vep] Done → $OUTPUT  ($SIZE)"
echo "[vep] Variants annotated:  $N"
echo "[vep] With gnomAD_min hit: $GNOMAD_HITS"
echo "════════════════════════════════════════════════════"
