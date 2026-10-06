#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────
#  Refresh ClinVar: the VEP custom VCF AND the PS1/PM5 amino-acid tables, from ONE release.
#
#  Usage:  update_clinvar.sh [clinvar.vcf.gz]
#          (no argument -> download the current NCBI weekly GRCh38 VCF + md5 check)
#
#  Writes, under $VEP_REFS/clinvar/ (dated, so earlier deliverables stay reproducible):
#    clinvar.chr.<fileDate>.vcf.gz(.tbi)   chr-named, primary contigs only
#    aa-<fileDate>/clinvar.MANE_missense.{PLP,BLB}.tsv + RELEASE
#  and repoints two symlinks at the new release:
#    clinvar.chr.vcf.gz(.tbi) -> clinvar.chr.<fileDate>.vcf.gz   ($CLINVAR_VCF default)
#    aa -> aa-<fileDate>                                          (set CLINVAR_AA_DIR to it)
#
#  The amino-acid tables are derived from THE SAME VCF that VEP annotates with, so PS1/PM5
#  and the per-variant ClinVar columns can no longer come from different releases.
#  Missense SNVs classified P/LP or B/LB are mapped to the MANE Select protein by VEP
#  (offline cache), giving gene / residue / ref AA / alt AA. Conflicting, VUS, somatic-only
#  ('.') and non-germline terms are left out, as in the previous variant_summary build.
#
#  Column layout matches filtering_r.pl load_clinvar_aa() (1-based):
#    3-6 Chr,PositionVCF,Ref,Alt  8 GeneSymbol  9 ClinicalSignificance  11 ReviewStatus
#    23 BB_AApos  24 BB_RefAA  28 AltAA  (other columns kept for layout, '.' when unused)
# ──────────────────────────────────────────────────────────────────────────
set -euo pipefail
. "$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)/site.sh"
export PERL5LIB="$PERL5LIB_EXTRA:${PERL5LIB:-}"
[ -e "$HTSLIB_SO" ] && export LD_PRELOAD="$HTSLIB_SO${LD_PRELOAD:+:$LD_PRELOAD}"

DEST="$VEP_REFS/clinvar"
URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
[[ -x "$VEP" ]] || { echo "ERROR: VEP not found at $VEP" >&2; exit 1; }
mkdir -p "$DEST"
TMP=$(mktemp -d "$DEST/.update.XXXXXX"); trap 'rm -rf "$TMP"' EXIT

# ── 1. Obtain the release ──
if [[ $# -ge 1 ]]; then
    RAW="$1"; [[ -s "$RAW" ]] || { echo "ERROR: $RAW not found" >&2; exit 1; }
    [[ -s "$RAW.tbi" ]] || tabix -f -p vcf "$RAW"
else
    RAW="$TMP/clinvar.vcf.gz"
    curl -fsS -o "$RAW" "$URL" -o "$RAW.tbi" "$URL.tbi" -o "$TMP/md5" "$URL.md5"
    echo "$(awk '{print $1}' "$TMP/md5")  $RAW" | md5sum -c --quiet \
        || { echo "ERROR: ClinVar md5 mismatch" >&2; exit 1; }
fi
FDATE=$(bcftools view -h "$RAW" | sed -n 's/^##fileDate=//p' | head -1)
[[ "$FDATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || { echo "ERROR: no ##fileDate in $RAW" >&2; exit 1; }
echo "[clinvar] release fileDate=$FDATE"

# ── 2. chr-named custom VCF (primary contigs only) ──
CHRVCF="$DEST/clinvar.chr.$FDATE.vcf.gz"
MAP="$TMP/chr_map.txt"
for c in {1..22} X Y; do printf '%s\tchr%s\n' "$c" "$c"; done > "$MAP"; printf 'MT\tchrM\n' >> "$MAP"
bcftools view -r "$(cut -f1 "$MAP" | paste -sd,)" "$RAW" -Ou \
  | bcftools annotate --rename-chrs "$MAP" -Oz -o "$TMP/chr.vcf.gz"
tabix -f -p vcf "$TMP/chr.vcf.gz"
mv -f "$TMP/chr.vcf.gz" "$CHRVCF"; mv -f "$TMP/chr.vcf.gz.tbi" "$CHRVCF.tbi"

# ── 3. P/LP + B/LB missense SNVs -> VEP (MANE Select protein coordinates) ──
# Classification bucket from the first '|'-term of CLNSIG (as the old ';'-term prefix match):
#   PLP: Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic (+ ',_low_penetrance')
#   BLB: Benign, Likely_benign, Benign/Likely_benign
bcftools view -i 'INFO/CLNVC="single_nucleotide_variant" && INFO/MC~"missense_variant"' "$CHRVCF" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\n' \
  | awk -F'\t' -v OFS='\t' '{
        s=$6; sub(/\|.*/,"",s)
        if (s ~ /^(Pathogenic|Likely_pathogenic)/) b="PLP"
        else if (s ~ /^(Benign|Likely_benign)/)    b="BLB"
        else next
        print $1,$2,$3,$4,$5,b,$6,$7 }' > "$TMP/sel.tsv"
echo "[clinvar] P/LP + B/LB missense SNVs: $(wc -l < "$TMP/sel.tsv")"
{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk -F'\t' -v OFS='\t' '{print $1,$2,$5,$3,$4,".",".","."}' "$TMP/sel.tsv"; } > "$TMP/sel.vcf"

"$VEP" --input_file "$TMP/sel.vcf" --output_file "$TMP/vep.tsv" --tab --no_stats \
  --cache --offline --dir_cache "$VEP_DATA" --assembly GRCh38 --species homo_sapiens \
  --fork "${VEP_FORKS:-4}" --buffer_size 5000 --mane_select --symbol --force_overwrite \
  --fields "Uploaded_variation,SYMBOL,Feature,Consequence,Protein_position,Amino_acids,MANE_SELECT"

# ── 4. Join + split into the two tables ──
AA="$DEST/aa-$FDATE"; mkdir -p "$AA"
HDR="VARKEY	KEY	Chr	PositionVCF	ReferenceAlleleVCF	AlternateAlleleVCF	Name	GeneSymbol	ClinicalSignificance	LastEvaluated	ReviewStatus	NumberSubmitters	PhenotypeIDs	PhenotypeList	VariationID	BB_Gene	BB_Transcript	BB_Strand	BB_CodingPOS	BB_CodingAllele	BB_CodonPOS	BB_CodonSeq	BB_AApos	BB_RefAA	ALT_Coding	RefCodon	AltCodon	AltAA	Consequence	MANE_HGVSc	MANE_HGVSp_1	MANE_HGVSp_3"
for b in PLP BLB; do printf '%s\n' "$HDR" > "$TMP/$b.tsv"; done
awk -F'\t' -v OFS='\t' -v T="$TMP" '
  NR==FNR { sel[$5]=$0; next }                          # ClinVar VariationID -> record
  /^#/ { next }
  $7=="" || $7=="-" || $4 !~ /missense_variant/ { next } # MANE Select missense only
  {
    if (!($1 in sel)) next
    split(sel[$1], s, "\t"); split($6, aa, "/")
    if (aa[1]=="" || aa[2]=="" || aa[1]==aa[2] || $5 !~ /^[0-9]+$/) next
    k=s[1]":"s[2]":"s[3]":"s[4]
    if (k in done) next; done[k]=1                      # one MANE record per variant
    D="."
    print k, s[1]":"s[2]":"s[3], s[1], s[2], s[3], s[4], D, $2, s[7], D, s[8], D, D, D, s[5], \
          $2, $3, D, D, D, D, D, $5, aa[1], D, D, D, aa[2], "missense", D, D, D  >> (T "/" s[6] ".tsv")
  }' "$TMP/sel.tsv" "$TMP/vep.tsv"
mv -f "$TMP/PLP.tsv" "$AA/clinvar.MANE_missense.PLP.tsv"
mv -f "$TMP/BLB.tsv" "$AA/clinvar.MANE_missense.BLB.tsv"
printf 'ClinVar_fileDate=%s\nsource=%s\nbuilt=%s\n' "$FDATE" "${1:-$URL}" "$(date -u +%FT%TZ)" > "$AA/RELEASE"

# ── 5. Point the defaults at the new release ──
ln -sfn "clinvar.chr.$FDATE.vcf.gz"     "$DEST/clinvar.chr.vcf.gz"
ln -sfn "clinvar.chr.$FDATE.vcf.gz.tbi" "$DEST/clinvar.chr.vcf.gz.tbi"
ln -sfn "aa-$FDATE"                      "$DEST/aa"
echo "[clinvar] PLP rows: $(($(wc -l < "$AA/clinvar.MANE_missense.PLP.tsv")-1))  BLB rows: $(($(wc -l < "$AA/clinvar.MANE_missense.BLB.tsv")-1))"
echo "[clinvar] done: $DEST/clinvar.chr.vcf.gz -> $FDATE ; set CLINVAR_AA_DIR=$DEST/aa"
