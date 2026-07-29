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
set -uo pipefail

REPO="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"
cd "$REPO"
for t in perl bgzip; do command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH"; exit 1; }; done

fails=0
ok()  { echo "  PASS  $*"; }
bad() { echo "  FAIL  $*"; fails=$((fails+1)); }

# ─────────────── 1 + 2: the built-in self-tests ───────────────
echo "== 1. family auto-discovery =="
if perl filtering_r.pl --selftest 2>&1 | grep -q 'all 9 passed'; then ok "naming self-test: 9/9"
else bad "naming self-test"; perl filtering_r.pl --selftest 2>&1 | tail -5; fi

echo "== 2. cohort recurrent-artifact filter =="
if perl filtering_r.pl --selftest-cohort 2>&1 | grep -q 'all 11 passed'; then ok "cohort self-test: 11/11"
else bad "cohort self-test"; perl filtering_r.pl --selftest-cohort 2>&1 | tail -5; fi

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
csq() { # <gene> <consequence> <cadd> <am_score> <ac> <an>
  printf '%s|1|%s|%s|%s|c.100G>A|p.Gly34Ser|100|G/S|34||||%s|likely_pathogenic|%s||||%s|%s|0|0|PASS|||' \
         "$1" "$MANE_TX" "$MANE_TX" "$2" "$3" "$4" "$5" "$6"
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
    [ "$nrow" -eq 1 ] && ok "exactly 1 candidate row (common + synonymous also dropped)" \
                      || bad "expected 1 row, got $nrow (genes: $(awk -F'\t' 'NR>1{print $6}' "$OUT" | tr '\n' ' '))"
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

echo
if [ "$fails" -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "$fails TEST(S) FAILED"; exit 1; fi
