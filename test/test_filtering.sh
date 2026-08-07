#!/usr/bin/env bash
#
# test_filtering.sh — regression test for the filtering logic. Synthetic data only:
# no patient data, no VEP, no cloud, no network, no GPU. Runs in a temp dir in ~5s.
# Requires: perl, bgzip (htslib). Run it after changing filtering_r.pl.
#
# Covers:
#   1. Family auto-discovery          (filtering_r.pl --selftest,        9 assertions)
#   2. Cohort recurrent-artifact gate (filtering_r.pl --selftest-cohort, 11 assertions)
#   3. Reference data integrity       — the tracked panel/whitelist/constraint files
#                                       still parse with the expected shape
#   4. End-to-end gating              — a synthetic annotated VCF through both passes:
#                                       the one real candidate is kept, and a common
#                                       variant, a synonymous one, a non-panel gene and a
#                                       solitary AR-gene het (carrier) are all dropped.
#   5. Consult mode + BP7             — --lookup output is namespaced, BP7 needs a real score
#   6. Batch-level table              — batch.<panel>.candidatos agrees with the per-proband one
#   7. Genotype + criterion edges     — haploid (hemizygous) GT, plain-XL MOI, the PP5 review-star
#                                       gate, per-gene row collapse, clinvar_conflict flagging
set -uo pipefail

REPO="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"
cd "$REPO"
for t in perl bgzip; do command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH"; exit 1; }; done

fails=0
ok()  { echo "  PASS  $*"; }
bad() { echo "  FAIL  $*"; fails=$((fails+1)); }

# ─────────────── 1 + 2: the built-in self-tests ───────────────
# Match "all N passed" without pinning N, so adding a case to a self-test does not
# fail this wrapper; the count is reported, and any [FAIL] line trips it.
selftest() { # <flag> <label>
    local out; out=$(perl filtering_r.pl "$1" 2>&1)
    local n; n=$(printf '%s' "$out" | sed -n 's/.*all \([0-9]\{1,\}\) passed.*/\1/p' | tail -1)
    if [ -n "$n" ] && ! printf '%s' "$out" | grep -q '\[FAIL\]'; then ok "$2: $n/$n"
    else bad "$2"; printf '%s\n' "$out" | tail -5; fi
}

echo "== 1. family auto-discovery =="
selftest --selftest        "naming self-test"

echo "== 2. cohort recurrent-artifact filter =="
selftest --selftest-cohort "cohort self-test"

# ─────────────── 3: reference data integrity ───────────────
# These files ARE the clinical behaviour — a truncated download silently changes results.
echo "== 3. tracked reference data =="
n=$(grep -vcE '^#|^$' g4e-2026.txt);            [ "$n" -gt 500 ]  && ok "panel: $n genes"            || bad "panel too small ($n)"
n=$(awk -F'\t' '!/^#/ && NF==4' g4e-2026.txt | wc -l)
[ "$n" -gt 500 ] && ok "panel: $n rows have the 4 expected columns" || bad "panel column shape ($n 4-col rows)"
n=$(grep -c . typevar.txt);                     [ "$n" -ge 8 ]    && ok "consequence whitelist: $n"  || bad "typevar too small ($n)"
n=$(grep -c '^ENST' mane-plus-clinical-names.txt); [ "$n" -gt 15000 ] && ok "MANE transcripts: $n"   || bad "MANE list too small ($n)"
n=$(grep -vcE '^#|^$' acmg_sf_v3.2.txt);        [ "$n" -eq 81 ]   && ok "ACMG SF genes: $n"          || bad "ACMG SF list wrong size ($n, expected 81)"
n=$(grep -c . gnomad-mis-constraint.txt);       [ "$n" -gt 15000 ] && ok "missense constraint: $n"   || bad "constraint table too small ($n)"

# ─────────────── 4: end-to-end gating on synthetic data ───────────────
echo "== 4. end-to-end (synthetic annotated VCF) =="
TD="$(mktemp -d)"; trap 'rm -rf "$TD"' EXIT
for f in filtering_r.pl parse_pangolin.pl g4e-2026.txt typevar.txt \
         mane-plus-clinical-names.txt acmg_sf_v3.2.txt gnomad-mis-constraint.txt; do
    ln -sf "$REPO/$f" "$TD/$f"
done
MANE_TX=$(grep -m1 '^ENST' "$REPO/mane-plus-clinical-names.txt" | cut -f1)
# A DOMINANT panel gene for the keep-case: a solitary het in a pure AR gene is
# deliberately dropped as a carrier, so an AR gene would test the wrong thing.
PANEL_GENE=$(awk -F'\t' '!/^#/ && $3=="AD" {print $1; exit}' "$REPO/g4e-2026.txt")
AR_GENE=$(awk -F'\t' '!/^#/ && $3=="AR" {print $1; exit}' "$REPO/g4e-2026.txt")
[ -n "$PANEL_GENE" ] && [ -n "$AR_GENE" ] || { echo "  FAIL  could not pick AD/AR genes from the panel"; exit 1; }

CSQ='SYMBOL|STRAND|Feature|MANE_SELECT|Consequence|HGVSc|HGVSp|cDNA_position|Amino_acids|Protein_position|REVEL|EVE_CLASS|EVE_SCORE|CADD_PHRED|am_class|am_pathogenicity|LoF|LoF_filter|LoF_flags|gnomADmin_AC_joint|gnomADmin_AN_joint|gnomADmin_AF_joint|gnomADmin_nhomalt_joint|gnomADmin_FILTER|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN'
# one CSQ record: gene, consequence, CADD, AM score, gnomAD AC/AN
csq() { # <gene> <consequence> <cadd> <am_score> <ac> <an> [clnsig] [clnrevstat]
  printf '%s|1|%s|%s|%s|c.100G>A|p.Gly34Ser|100|G/S|34||||%s|likely_pathogenic|%s||||%s|%s|0|0|PASS|%s|%s|' \
         "$1" "$MANE_TX" "$MANE_TX" "$2" "$3" "$4" "$5" "$6" "${7:-}" "${8:-}"
}
STAR1='criteria_provided,_single_submitter'
# A dual-inheritance gene (MOI lists both AD and AR): its solitary hets pass through
# as dominant candidates, which is what the recessive-flag row-consistency check needs.
DUAL_GENE=$(awk -F'\t' '!/^#/ && $3 ~ /AD/ && $3 ~ /AR/ {print $1; exit}' "$REPO/g4e-2026.txt")
[ -n "$DUAL_GENE" ] || { echo "  FAIL  no dual-inheritance gene in the panel"; exit 1; }
{
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr2,length=250000000>'
  echo '##FILTER=<ID=PASS,Description="p">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">'
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">'
  echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">'
  echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">'
  echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: $CSQ\">"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTFAM-P\n'
  # KEEP: rare, panel gene, MANE, missense, CADD 30, AlphaMissense 0.99, absent from gnomAD
  printf 'chr2\t1000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" missense_variant 30 0.99 0 0)"
  # DROP: same but COMMON in gnomAD (5%)
  printf 'chr2\t2000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" missense_variant 30 0.99 5000 100000)"
  # DROP: synonymous (not on the consequence whitelist)
  printf 'chr2\t3000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" synonymous_variant 30 0.99 0 0)"
  # DROP: gene not on the panel (and not an ACMG SF gene)
  printf 'chr2\t4000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq ZZZNOTAGENE missense_variant 30 0.99 0 0)"
  # DROP: solitary het in a pure AR gene = carrier state, not an explanation (default drop)
  printf 'chr2\t5000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$AR_GENE" missense_variant 30 0.99 0 0)"
  # KEEP: ClinVar P/LP (1 star) at 5% in gnomAD — an established classification outranks
  # the rarity ceiling, so founder alleles are not lost. No in-silico support at all.
  printf 'chr2\t6000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" missense_variant 5 0.10 5000 100000 Pathogenic "$STAR1")"
  # KEEP: ClinVar P/LP (1 star) on a consequence that is NOT whitelisted (intronic)
  printf 'chr2\t7000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" intron_variant 5 0.10 0 0 Likely_pathogenic "$STAR1")"
  # KEEP via LoF: 2 bp deletion, used to check the END coordinate spans the REF allele
  printf 'chr2\t8000\t.\tCAA\tC\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" frameshift_variant 5 0.10 0 0)"
  # KEEP: homozygous + heterozygous variant in the SAME dual-inheritance gene. The gene
  # is flagged HOM; the het row must not inherit that label.
  printf 'chr2\t9000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t1/1:0,40:40:99\n' "$(csq "$DUAL_GENE" missense_variant 30 0.99 0 0)"
  printf 'chr2\t9500\t.\tC\tT\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$DUAL_GENE" missense_variant 30 0.99 0 0)"
  # KEEP: homozygous variant in a purely DOMINANT gene. Recessive reasoning does not
  # apply to such genes, so it must carry no recessive flag despite being hom.
  printf 'chr2\t10000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t1/1:0,40:40:99\n' "$(csq "$PANEL_GENE" missense_variant 30 0.99 0 0)"
} | bgzip -c > "$TD/TESTFAM-P.germline.vep.vcf.gz"

# Pass 1 emits <proband>.<panel>.pangolin_input.csv and stops. Supply an empty score map
# (same panel-namespaced prefix) so pass 2 runs without Pangolin/GPU — variants that need a
# splice score simply get none, which is the documented graceful behaviour.
( cd "$TD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >pass1.log 2>&1 )
IN=$(ls "$TD"/TESTFAM-P.*.pangolin_input.csv 2>/dev/null | head -1)
[ -n "$IN" ] && ok "pass 1 emitted $(basename "$IN")" || bad "pass 1 emitted no pangolin input"
: > "${IN%.pangolin_input.csv}.pangolin.tsv"
( cd "$TD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >pass2.log 2>&1 )

OUT=$(ls "$TD"/TESTFAM-P.*.candidatos 2>/dev/null | head -1)
if [ -z "$OUT" ]; then
    bad "no .candidatos produced"; tail -15 "$TD/pass2.log" | sed 's/^/      /'
else
    ok "produced $(basename "$OUT")"
    kept=$(awk -F'\t' 'NR>1{print $2}' "$OUT" | sort -u | tr '\n' ' ')
    nrow=$(awk 'NR>1' "$OUT" | wc -l)
    grep -q "$PANEL_GENE" "$OUT" && ok "rare damaging panel variant KEPT ($PANEL_GENE)" \
                                 || bad "rare damaging panel variant was dropped"
    grep -q 'ZZZNOTAGENE' "$OUT" && bad "non-panel gene leaked into output" \
                                 || ok "non-panel gene dropped"
    grep -q "	$AR_GENE	" "$OUT" && bad "solitary AR-gene het (carrier) leaked into output" \
                                  || ok "solitary AR-gene het dropped (carrier rule)"
    row() { awk -F'\t' -v p="$1" 'NR>1 && $2==p' "$OUT"; }          # a data row by POS
    col() { awk -F'\t' -v p="$1" -v c="$2" 'NR>1 && $2==p{print $c}' "$OUT"; }
    # Resolve a column by HEADER NAME, the way filtering_r.pl resolves CSQ fields.
    # Hard-coded indices silently pointed at a neighbouring column when the table
    # changed shape (recessive_flag + qc_flag became one `flags` column).
    ncol() { awk -F'\t' -v n="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==n){print i; exit}}' "$OUT"; }
    # One ';'-separated component of the consolidated flags field, by POS.
    hasflag() { awk -F'\t' -v p="$1" -v c="$2" -v f="$3" \
                    'NR>1 && $2==p{n=split($c,a,";"); for(i=1;i<=n;i++) if(a[i]==f){print "yes"; exit}}' "$OUT"; }

    # An established ClinVar classification outranks the rarity ceiling and the
    # consequence whitelist — both gates otherwise fire before any evidence is read.
    [ -n "$(row 6000)" ] && ok "ClinVar P/LP kept at 5% gnomAD AF (founder allele preserved)" \
                         || bad "ClinVar P/LP variant dropped by the rarity gate"
    [ -n "$(row 7000)" ] && ok "ClinVar P/LP kept on a non-whitelisted consequence (intronic)" \
                         || bad "ClinVar P/LP variant dropped by the consequence whitelist"

    e=$(col 8000 3)
    [ "$e" = "8002" ] && ok "indel END spans the REF allele (8000 CAA/C -> 8002)" \
                      || bad "indel END wrong: got '${e:-none}', expected 8002"

    FC=$(ncol flags); ZC=$(ncol zygosity)
    [ -n "$FC" ] && ok "consolidated 'flags' column present in the header" \
                 || bad "no 'flags' column in the output header"

    if [ -n "$(row 9500)" ] && [ -n "$(row 9000)" ]; then
        [ -z "$(hasflag 9500 "$FC" HOM)" ] \
            && ok "het row not labelled HOM by a hom variant elsewhere in the gene" \
            || bad "het row carries HOM in flags while zygosity=$(col 9500 "$ZC")"
        [ -n "$(hasflag 9000 "$FC" HOM)" ] \
            && ok "hom row still labelled HOM" \
            || bad "hom row lost its HOM label (flags='$(col 9000 "$FC")')"
    else
        bad "dual-inheritance gene rows missing from output"
    fi

    dh=$(col 10000 "$FC")
    [ -n "$(row 10000)" ] && [ -z "$(hasflag 10000 "$FC" HOM)" ] \
        && ok "hom row in a purely dominant gene carries no recessive flag" \
        || bad "dominant-gene hom row got flags='${dh:-?}' (expected no HOM)"

    [ "$nrow" -eq 7 ] && ok "exactly 7 candidate rows (common, synonymous, off-panel, carrier dropped)" \
                      || bad "expected 7 rows, got $nrow (POS: $(awk -F'\t' 'NR>1{print $2}' "$OUT" | tr '\n' ' '))"
fi

# ─────────────── 5: consult mode + BP7 splice evidence ───────────────
echo "== 5. consult mode (--lookup) and BP7 =="
if [ -n "${OUT:-}" ]; then
    TSV="${IN%.pangolin_input.csv}.pangolin.tsv"
    before=$(md5sum "$OUT" | cut -d' ' -f1)
    ( cd "$TD" && CLINVAR_AA_DIR= REF_FASTA= \
        perl filtering_r.pl --lookup TESTFAM-P.germline.vep.vcf.gz >lookup.log 2>&1 )
    LK=$(ls "$TD"/Lookup.TESTFAM-P.*.candidatos 2>/dev/null | head -1)
    [ -n "$LK" ] && ok "--lookup wrote $(basename "$LK")" \
                 || { bad "--lookup produced no Lookup.* output"; tail -8 "$TD/lookup.log" | sed 's/^/      /'; }
    # The consult must never clobber the cohort table: --lookup derives its tag
    # from the input VCF basename, which IS the proband name.
    [ "$(md5sum "$OUT" | cut -d' ' -f1)" = "$before" ] \
        && ok "cohort candidatos table left untouched by --lookup" \
        || bad "--lookup overwrote $(basename "$OUT")"

    # BP7 asserts ABSENCE of splice impact. With no score map it must stay silent;
    # with a real sub-threshold score it must fire. Both directions are checked so
    # the guard cannot be satisfied by simply disabling the criterion.
    if [ -n "$LK" ]; then
        grep -q 'BP7' "$LK" && bad "BP7 asserted with no Pangolin score" \
                            || ok "BP7 withheld when splicing was never scored"
    fi
    printf 'chr2-3000-G-A\t0.05\n' > "$TSV"
    ( cd "$TD" && CLINVAR_AA_DIR= REF_FASTA= \
        perl filtering_r.pl --lookup TESTFAM-P.germline.vep.vcf.gz >lookup2.log 2>&1 )
    LK2=$(ls "$TD"/Lookup.TESTFAM-P.*.candidatos 2>/dev/null | head -1)
    grep -q 'BP7' "$LK2" && ok "BP7 fires on a scored synonymous variant (0.05 < 0.2)" \
                         || bad "BP7 did not fire despite a sub-threshold Pangolin score"
else
    bad "skipping consult tests — section 4 produced no candidatos"
fi

# ─────────────── 6: batch-level table ───────────────
echo "== 6. batch-level table =="
BATCH=$(ls "$TD"/batch.*.candidatos 2>/dev/null | head -1)
if [ -n "$BATCH" ]; then
    ok "wrote $(basename "$BATCH")"
    [ "$(head -1 "$BATCH" | cut -f1)" = "sample" ] \
        && ok "batch table's first column is 'sample'" \
        || bad "batch table's first column is '$(head -1 "$BATCH" | cut -f1)', expected 'sample'"
    # Same rows as the per-proband table, plus the sample prefix — the two must agree,
    # otherwise the batch view and the per-case view tell a curator different stories.
    nb=$(awk 'NR>1' "$BATCH" | wc -l)
    [ "$nb" -eq "${nrow:-0}" ] && ok "batch rows match the per-proband table ($nb)" \
                              || bad "batch has $nb rows, per-proband table has ${nrow:-0}"
    [ "$(awk -F'\t' 'NR>1{print $1}' "$BATCH" | sort -u)" = "TESTFAM-P" ] \
        && ok "batch rows carry the proband id" \
        || bad "unexpected sample ids in batch table"
    [ "$(head -1 "$BATCH" | cut -f2-)" = "$(head -1 "$OUT")" ] \
        && ok "batch columns 2..N are identical to the per-proband header" \
        || bad "batch header diverges from the per-proband header"
else
    bad "no batch.*.candidatos produced"
fi

# ─────────────── 7: haploid genotypes, X-linked genes, PP5 star gate ───────────────
# Everything here is a genotype/criterion case that a WGS/DRAGEN run hits routinely but
# the synthetic autosomal fixture above cannot reach.
echo "== 7. haploid GT, X-linked MOI, PP5 star gate, per-gene collapse =="
XD="$(mktemp -d)"; trap 'rm -rf "$TD" "$XD"' EXIT
for f in filtering_r.pl parse_pangolin.pl g4e-2026.txt typevar.txt \
         mane-plus-clinical-names.txt acmg_sf_v3.2.txt gnomad-mis-constraint.txt; do
    ln -sf "$REPO/$f" "$XD/$f"
done
# A plain "XL" gene — the g4e vocabulary for X-linked genes with no XLD/XLR split.
XL_GENE=$(awk -F'\t' '!/^#/ && $3=="XL" {print $1; exit}' "$REPO/g4e-2026.txt")
TX2=$(grep -m2 '^ENST' "$REPO/mane-plus-clinical-names.txt" | tail -1 | cut -f1)
GENE_B=$(awk -F'\t' '!/^#/ && $3=="AD" {print $1}' "$REPO/g4e-2026.txt" | sed -n 2p)
[ -n "$XL_GENE" ] && [ -n "$GENE_B" ] || { echo "  FAIL  could not pick XL/second AD gene"; exit 1; }
NOSTAR='no_assertion_criteria_provided'
csq2() { # like csq() but with an explicit transcript as $9
  printf '%s|1|%s|%s|%s|c.100G>A|p.Gly34Ser|100|G/S|34||||%s|likely_pathogenic|%s||||%s|%s|0|0|PASS|%s|%s|' \
         "$1" "${9:-$MANE_TX}" "${9:-$MANE_TX}" "$2" "$3" "$4" "$5" "$6" "${7:-}" "${8:-}"
}
xhdr() {
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr2,length=250000000>'
  echo '##contig=<ID=chrX,length=156040895>'
  echo '##FILTER=<ID=PASS,Description="p">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">'
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">'
  echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">'
  echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">'
  echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: $CSQ\">"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1"
}
{
  xhdr XFAM-P
  # HEMIZYGOUS male call in an X-linked gene: DRAGEN emits a single-allele GT.
  # No in-silico support at all — it must survive on the AR_hem rescue alone.
  printf 'chrX\t1000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t1:0,40:40:99\n' "$(csq "$XL_GENE" missense_variant 5 0.10 0 0)"
  # ClinVar Pathogenic with ZERO review stars: kept by the ClinVar rescue arm, but
  # must NOT earn PP5 (which now requires >=1 star, like every other ClinVar consumer).
  printf 'chr2\t2000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" missense_variant 5 0.10 0 0 Pathogenic "$NOSTAR")"
  # One position annotated onto TWO overlapping MANE gene models. Both are panel genes
  # and both pass the gates, so both rows must survive the one-row-per-variant collapse.
  printf 'chr2\t3000\t.\tG\tA\t500\tPASS\tCSQ=%s,%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
         "$(csq2 "$PANEL_GENE" missense_variant 30 0.99 0 0 '' '' "$MANE_TX")" \
         "$(csq2 "$GENE_B"     missense_variant 30 0.99 0 0 '' '' "$TX2")"
  # Truncating + gnomAD-absent (PVS1,PM2 -> Likely_pathogenic) but ClinVar says BENIGN
  # with a review star (BP6). The combining rules cannot express this, so the class
  # stays LP and the row must be flagged clinvar_conflict instead.
  printf 'chr2\t4000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq "$PANEL_GENE" stop_gained 35 0.99 0 0 Benign "$STAR1")"
} | bgzip -c > "$XD/XFAM-P.germline.vep.vcf.gz"
# Father: hemizygous for the same chrX variant. He can only be seen as a carrier if
# haploid GTs parse — otherwise the son's variant is called de novo.
{ xhdr XFAM-F
  printf 'chrX\t1000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t1:0,40:40:99\n' "$(csq "$XL_GENE" missense_variant 5 0.10 0 0)"
} | bgzip -c > "$XD/XFAM-F.germline.vep.vcf.gz"

( cd "$XD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >xpass1.log 2>&1 )
XIN=$(ls "$XD"/XFAM-P.*.pangolin_input.csv 2>/dev/null | head -1)
[ -n "$XIN" ] && : > "${XIN%.pangolin_input.csv}.pangolin.tsv"
( cd "$XD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >xpass2.log 2>&1 )
XOUT=$(ls "$XD"/XFAM-P.*.candidatos 2>/dev/null | head -1)

if [ -z "$XOUT" ]; then
    bad "section 7 produced no candidatos"; tail -15 "$XD/xpass2.log" | sed 's/^/      /'
else
    xcol() { awk -F'\t' -v n="$1" -v p="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $2==p{print $h[n]}' "$XOUT"; }
    z=$(xcol zygosity 1000)
    [ "$z" = "hem" ] && ok "haploid GT=1 reported as zygosity 'hem'" \
                     || bad "haploid GT=1 gave zygosity '${z:-empty}' (expected hem)"
    k=$(xcol kept_by 1000)
    case "$k" in *AR_hem*) ok "hemizygous X variant rescued by AR_hem (XL gene treated as recessive-capable)";;
                 *) bad "hemizygous X variant not rescued: kept_by='${k:-none}'";; esac
    f=$(xcol flags 1000)
    case ";$f;" in *";HEM;"*) ok "hemizygous row flagged HEM";;
                   *) bad "hemizygous row flags='${f:-empty}' (expected HEM)";; esac
    inh=$(xcol inheritance 1000)
    [ "$inh" = "IF" ] && ok "hemizygous father recognised as carrier (inheritance=IF, not de novo)" \
                      || bad "inheritance='${inh:-empty}' (expected IF; a hemizygous parent was invisible)"

    crit=$(xcol acmg_criteria 2000)
    [ -n "$(awk -F'\t' 'NR>1 && $2==2000' "$XOUT")" ] \
        && ok "0-star ClinVar P/LP still kept by the ClinVar rescue arm" \
        || bad "0-star ClinVar P/LP variant was dropped"
    case ",$crit," in *,PP5,*) bad "PP5 fired on a 0-star ClinVar submission (criteria=$crit)";;
                      *) ok "PP5 withheld on a 0-star ClinVar submission";; esac

    ng=$(awk -F'\t' 'NR>1 && $2==3000' "$XOUT" | wc -l)
    [ "$ng" -eq 2 ] && ok "overlapping gene models both survive the collapse (2 rows at one position)" \
                    || bad "collapse kept $ng row(s) at the two-gene position, expected 2"

    cls=$(xcol acmg_class 4000); ccr=$(xcol acmg_criteria 4000); cfl=$(xcol flags 4000)
    case "$ccr" in *BP6*) ok "BP6 fired on the starred ClinVar-Benign truncating variant";;
                    *) bad "BP6 missing on a 1-star ClinVar Benign variant (criteria=$ccr)";; esac
    case "$cls" in *athogenic*) ok "contradictory row still auto-classed $cls (class left alone by design)";;
                   *) bad "expected a pathogenic-leaning class, got '$cls'";; esac
    case ";$cfl;" in *";clinvar_conflict;"*) ok "contradiction surfaced as flags=clinvar_conflict";;
                     *) bad "no clinvar_conflict flag on $cls with criteria=$ccr (flags='${cfl:-empty}')";; esac
fi

# ─────────── 8: splice discovery, ACMG-SF scoring, PM2 strength, PS2 scope ───────────
echo "== 8. splice discovery, SF scoring, PM2 strength, PS2 scope =="
SD="$(mktemp -d)"; trap 'rm -rf "$TD" "$XD" "$SD"' EXIT
for f in filtering_r.pl parse_pangolin.pl g4e-2026.txt typevar.txt \
         mane-plus-clinical-names.txt acmg_sf_v3.2.txt gnomad-mis-constraint.txt; do
    ln -sf "$REPO/$f" "$SD/$f"
done
SF_GENE=$(awk -F'\t' '!/^#/ && NF{print $1}' "$REPO/acmg_sf_v3.2.txt" | grep -vx TTN | head -1)
# csq with an explicit HGVSc, so an intron offset (c.100+50) can be expressed.
csq3() { # <gene> <consequence> <cadd> <am> <ac> <an> <hgvsc>
  printf '%s|1|%s|%s|%s|%s|p.Gly34Ser|100|G/S|34||||%s|likely_pathogenic|%s||||%s|%s|0|0|PASS|||' \
         "$1" "$MANE_TX" "$MANE_TX" "$2" "$7" "$3" "$4" "$5" "$6"
}
{
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr2,length=250000000>'
  echo '##FILTER=<ID=PASS,Description="p">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">'
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">'
  echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">'
  echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">'
  echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: $CSQ\">"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPFAM-P\n'
  # PROBE, in window (c.100+50): deep-intronic, gnomAD-covered and rare. Not whitelisted,
  # so ONLY a positive Pangolin score can rescue it. Scored high below -> must appear.
  printf 'chr2\t1000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" intron_variant 5 0.10 1 200000 'c.100+50G>A')"
  # PROBE, in window, but Pangolin says no -> must NOT appear.
  printf 'chr2\t1100\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" intron_variant 5 0.10 1 200000 'c.200+60G>A')"
  # PROBE out of window (c.100+9000): beyond $INTRON_MAX_DIST -> never even scored.
  printf 'chr2\t1200\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" intron_variant 5 0.10 1 200000 'c.100+9000G>A')"
  # PROBE with NO gnomAD coverage (AN=0): frequency unknowable and PM2 would be
  # manufactured, so it must not be probed at all.
  printf 'chr2\t1300\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" intron_variant 5 0.10 0 0 'c.300+40G>A')"
  # Synonymous probe, gnomAD-covered: exonic splice-altering synonymous. Scored high.
  printf 'chr2\t1400\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" synonymous_variant 5 0.10 1 200000 'c.400G>A')"
  # ACMG-SF gene, whitelisted consequence -> must reach the Pangolin input (#6).
  printf 'chr2\t2000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$SF_GENE" missense_variant 30 0.99 1 200000 'c.500G>A')"
  # Truncating + gnomAD-absent, with NO computational support (empty AlphaMissense/REVEL)
  # so the class rests on PVS1 + PM2 alone — that is what makes it a PM2-strength probe.
  printf 'chr2\t3000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' "$(csq3 "$PANEL_GENE" stop_gained 5 '' 0 0 'c.600G>A')"
} | bgzip -c > "$SD/SPFAM-P.germline.vep.vcf.gz"

( cd "$SD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >spass1.log 2>&1 )
SIN=$(ls "$SD"/SPFAM-P.*.pangolin_input.csv 2>/dev/null | head -1)
if [ -z "$SIN" ]; then
    bad "section 8: pass 1 emitted no pangolin input"; tail -10 "$SD/spass1.log" | sed 's/^/      /'
else
    inwin=$(awk -F, 'NR>1 && $2==1000' "$SIN" | wc -l)
    [ "$inwin" -eq 1 ] && ok "intronic probe inside the window entered the Pangolin set" \
                       || bad "in-window intronic probe not scored"
    out=$(awk -F, 'NR>1 && $2==1200' "$SIN" | wc -l)
    [ "$out" -eq 0 ] && ok "intronic variant beyond \$INTRON_MAX_DIST not probed" \
                     || bad "out-of-window intronic variant was probed"
    nocov=$(awk -F, 'NR>1 && $2==1300' "$SIN" | wc -l)
    [ "$nocov" -eq 0 ] && ok "gnomAD-uncovered intronic variant not probed (would manufacture PM2)" \
                       || bad "uncovered intronic variant was probed"
    syn=$(awk -F, 'NR>1 && $2==1400' "$SIN" | wc -l)
    [ "$syn" -eq 1 ] && ok "synonymous probe entered the Pangolin set" \
                     || bad "synonymous probe not scored"
    sfin=$(awk -F, 'NR>1 && $2==2000' "$SIN" | wc -l)
    [ "$sfin" -eq 1 ] && ok "ACMG-SF gene variant reaches the Pangolin input (#6)" \
                      || bad "ACMG-SF variant still absent from the Pangolin input"

    # Score: rescue 1000 and 1400, reject 1100.
    STSV="${SIN%.pangolin_input.csv}.pangolin.tsv"
    { printf 'chr2-1000-G-A\t0.90\n'; printf 'chr2-1100-G-A\t0.05\n'; printf 'chr2-1400-G-A\t0.80\n'; } > "$STSV"
    ( cd "$SD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >spass2.log 2>&1 )
    SOUT=$(ls "$SD"/SPFAM-P.*.candidatos 2>/dev/null | head -1)
    if [ -z "$SOUT" ]; then
        bad "section 8: no candidatos"; tail -10 "$SD/spass2.log" | sed 's/^/      /'
    else
        scol() { awk -F'\t' -v n="$1" -v p="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $2==p{print $h[n]}' "$SOUT"; }
        [ -n "$(awk -F'\t' 'NR>1 && $2==1000' "$SOUT")" ] \
            && ok "DEEP-INTRONIC variant DISCOVERED via Pangolin (kept_by=$(scol kept_by 1000))" \
            || bad "high-scoring deep-intronic probe was not rescued"
        [ -z "$(awk -F'\t' 'NR>1 && $2==1100' "$SOUT")" ] \
            && ok "low-scoring probe stays out of the table" \
            || bad "probe with Pangolin 0.05 leaked into the table"
        [ -n "$(awk -F'\t' 'NR>1 && $2==1400' "$SOUT")" ] \
            && ok "splice-altering SYNONYMOUS variant discovered" \
            || bad "high-scoring synonymous probe was not rescued"

        # PM2 strength: PVS1 + PM2 alone. At Supporting (ClinGen SVI) 2015 has no
        # "PVS1 + 1 supporting" rule, so this is VUS; at Moderate it would be LP.
        pm2c=$(scol acmg_criteria 3000); pm2k=$(scol acmg_class 3000)
        case ",$pm2c," in *,PM2_Supporting,*) ok "PM2 recorded as PM2_Supporting (ClinGen SVI)";;
                          *,PM2,*) ok "PM2 recorded at Moderate (ACMG 2015 reading)";;
                          *) bad "no PM2 on a gnomAD-absent variant (criteria=$pm2c)";; esac
        case "$pm2c" in
          *PM2_Supporting*) [ "$pm2k" = "VUS" ] \
              && ok "PVS1+PM2_Supporting combines to VUS (no such rule in ACMG 2015 Table 5)" \
              || bad "PVS1+PM2_Supporting gave '$pm2k', expected VUS" ;;
          *) [ "$pm2k" = "Likely_pathogenic" ] \
              && ok "PVS1+PM2(Moderate) combines to Likely_pathogenic" \
              || bad "PVS1+PM2 gave '$pm2k', expected Likely_pathogenic" ;;
        esac

        # An uncovered probe must not earn PM2: AN=0 there means "outside the resource",
        # not "unobserved". Re-run with --probe-uncovered so 1300 (AN=0) is probed, and
        # score it high enough to be rescued.
        printf 'chr2-1300-G-A\t0.95\n' >> "$STSV"
        ( cd "$SD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl --probe-uncovered >spass3.log 2>&1 )
        upm=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $2==1300{print $h["acmg_criteria"]}' "$SOUT")
        if [ -n "$(awk -F'\t' 'NR>1 && $2==1300' "$SOUT")" ]; then
            ok "--probe-uncovered widens the probe set past the gnomAD footprint"
            case ",$upm," in *,PM2*) bad "PM2 awarded on an uncovered probe (criteria=$upm)";;
                             *) ok "PM2 withheld on a gnomAD-uncovered probe row";; esac
        else
            bad "--probe-uncovered did not rescue the uncovered probe"
        fi

        # PS2 scope: this is a SINGLETON run, so inheritance is NA and PS2 must never fire.
        grep -q 'PS2' "$SOUT" && bad "PS2 fired in a singleton run (no parents present)" \
                              || ok "PS2 never fires without a trio (singleton run)"
    fi
fi

echo "== 9. PM1 from the PERv1 custom track =="
PD="$(mktemp -d)"; trap 'rm -rf "$TD" "$XD" "$SD" "$PD"' EXIT
for f in filtering_r.pl parse_pangolin.pl g4e-2026.txt typevar.txt \
         mane-plus-clinical-names.txt acmg_sf_v3.2.txt gnomad-mis-constraint.txt; do
    ln -sf "$REPO/$f" "$PD/$f"
done
# CSQ layout with the PER field appended, exactly as vep_annotate.sh emits it.
CSQP="$CSQ|PER"
# <gene> <consequence> <PER name field>
csq4() { printf '%s|1|%s|%s|%s|c.100G>A|p.Gly34Ser|100|G/S|34||||30|likely_pathogenic|0.99||||0|0|0|0|PASS||||%s' \
                "$1" "$MANE_TX" "$MANE_TX" "$2" "$3"; }
{
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr2,length=250000000>'
  echo '##FILTER=<ID=PASS,Description="p">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">'
  echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">'
  echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">'
  echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">'
  echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: $CSQP\">"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPMFAM-P\n'
  # 1000: missense in a Strong PER of THIS gene -> PM1_Strong.
  printf 'chr2\t1000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "$PANEL_GENE/PERv1_direct/PM1_Strong/PER1/p.30-40/OR=42.9/adjP=0.01")"
  # 2000: missense in a Moderate PER -> PM1 at Moderate.
  printf 'chr2\t2000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "$PANEL_GENE/PERv1_direct/PM1_Moderate/PER2/p.50-60/OR=6.1/adjP=0.02")"
  # 3000: SYNONYMOUS inside a PER. A PER is a missense-burden statement -> no PM1.
  printf 'chr2\t3000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" synonymous_variant "$PANEL_GENE/PERv1_direct/PM1_Strong/PER3/p.70-80/OR=99/adjP=0.001")"
  # 4000: missense overlapping a PER computed on a DIFFERENT gene (co-located MANE
  # models are real). The region is evidence about that gene, not this one -> no PM1.
  printf 'chr2\t4000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "OTHERGENE/PERv1_direct/PM1_Strong/PER1/p.90-99/OR=50/adjP=0.001")"
  # 5000: missense with NO PER overlap -> no PM1.
  printf 'chr2\t5000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "")"
  # 6000: missense in a PARALOG region only -> PM1 at Moderate (transfer costs a tier,
  # so it is capped even though this region's enrichment is high).
  printf 'chr2\t6000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "$PANEL_GENE/PERv1_paralog/PM1_Moderate/PER1/p.10-20/OR=99/adjP=1e-08")"
  # 7000: overlaps BOTH a direct-Strong and a paralog region. The direct one is
  # evidence about this gene rather than transferred onto it, so it must win.
  printf 'chr2\t7000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "$PANEL_GENE/PERv1_paralog/PM1_Moderate/PER9/p.1-9/OR=50/adjP=0.01&$PANEL_GENE/PERv1_direct/PM1_Strong/PER4/p.5-15/OR=30/adjP=0.002")"
  # 8000: paralog region belonging to ANOTHER gene -> no PM1 (gene guard).
  printf 'chr2\t8000\t.\tG\tA\t500\tPASS\tCSQ=%s\tGT:AD:DP:GQ\t0/1:20,20:40:99\n' \
    "$(csq4 "$PANEL_GENE" missense_variant "OTHERGENE/PERv1_paralog/PM1_Moderate/PER1/p.1-9/OR=50/adjP=0.01")"
} | bgzip -c > "$PD/PMFAM-P.germline.vep.vcf.gz"
# EMIT pass, then an empty Pangolin table so the FINAL pass can run.
( cd "$PD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >pm1.emit.log 2>&1 )
PIN=$(ls "$PD"/*.pangolin_input.csv 2>/dev/null | head -1)
[ -n "$PIN" ] && : > "${PIN%.pangolin_input.csv}.pangolin.tsv"
( cd "$PD" && CLINVAR_AA_DIR= REF_FASTA= perl filtering_r.pl >pm1.log 2>&1 )
POUT=$(ls "$PD"/PMFAM-P.*.candidatos 2>/dev/null | head -1)
if [ -z "$POUT" ]; then
    bad "section 9: no candidatos"; tail -10 "$PD/pm1.log" | sed 's/^/      /'
else
    pcol() { awk -F'\t' -v n="$1" -v p="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $2==p{print $h[n]}' "$POUT"; }
    case ",$(pcol acmg_criteria 1000)," in *,PM1_Strong,*) ok "PM1_Strong awarded in a high-enrichment PER";;
        *) bad "PM1_Strong not awarded (criteria=$(pcol acmg_criteria 1000))";; esac
    case "$(pcol flags 1000)" in *PM1:PERv1_direct/PER1/p.30-40*) ok "PM1 region + arm recorded in flags for audit";;
        *) bad "PM1 region detail missing from flags (flags=$(pcol flags 1000))";; esac
    case ",$(pcol acmg_criteria 2000)," in *,PM1,*) ok "PM1 awarded at Moderate in a low-enrichment PER";;
        *) bad "PM1 (Moderate) not awarded (criteria=$(pcol acmg_criteria 2000))";; esac
    case ",$(pcol acmg_criteria 3000)," in *,PM1*) bad "PM1 fired on a SYNONYMOUS variant";;
        *) ok "PM1 withheld on a synonymous variant inside a PER";; esac
    case ",$(pcol acmg_criteria 4000)," in *,PM1*) bad "PM1 fired from a PER belonging to another gene";;
        *) ok "PM1 withheld when the PER belongs to a different gene";; esac
    case ",$(pcol acmg_criteria 5000)," in *,PM1*) bad "PM1 fired without any PER overlap";;
        *) ok "PM1 withheld outside every PER";; esac
    case ",$(pcol acmg_criteria 6000)," in *,PM1,*) ok "paralog PER awards PM1 at Moderate";;
        *,PM1_Strong,*) bad "paralog PER awarded Strong (cap not applied)";;
        *) bad "paralog PER awarded no PM1 (criteria=$(pcol acmg_criteria 6000))";; esac
    case "$(pcol flags 6000)" in *PERv1_paralog*) ok "flags record which arm fired (paralog)";;
        *) bad "arm not recorded in flags (flags=$(pcol flags 6000))";; esac
    case ",$(pcol acmg_criteria 7000)," in *,PM1_Strong,*) ok "direct PER outranks an overlapping paralog PER";;
        *) bad "direct did not win over paralog (criteria=$(pcol acmg_criteria 7000))";; esac
    case "$(pcol flags 7000)" in *PERv1_direct*) ok "flags record the winning arm (direct)";;
        *) bad "wrong arm recorded (flags=$(pcol flags 7000))";; esac
    case ",$(pcol acmg_criteria 8000)," in *,PM1*) bad "PM1 fired from a paralog PER of another gene";;
        *) ok "PM1 withheld when the paralog PER belongs to a different gene";; esac
fi

echo
if [ "$fails" -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "$fails TEST(S) FAILED"; exit 1; fi
