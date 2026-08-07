#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

# Naming / family-discovery self-test (no data or reference files needed).
#   perl filtering_r.pl --selftest
run_naming_selftest() if grep { $_ eq '--selftest' } @ARGV;

#############################################################################
# filtering_r.pl  —  clinical candidate filtering for trio/duo germline VCFs
#                    (robust successor of filtering_new.pl / filtering_b3.pl)
#
# Pipeline features
# -----------------
#  * CSQ fields resolved BY NAME from each VCF's own CSQ header (no hard-coded
#    indices). Critical fields are asserted at startup (loud failure, never
#    silent).                                                          [#3]
#  * Reads *.germline.vep.vcf.gz directly; auto-discovers families by filename
#    using the role-suffix convention (see sample_role): <FAMILY>-P=proband,
#    <FAMILY>-M=mother, <FAMILY>-F=father. A name with no -P/-M/-F suffix is
#    ignored by discovery; if NO input has one, every sample is analyzed as a
#    singleton proband.
#  * Multiallelic input is assumed pre-split (vep_annotate.sh runs
#    `bcftools norm -m-any`). Any residual multiallelic row is skipped+counted.
#                                                                      [#2]
#  * Structural gates (AND): MANE transcript + whitelisted consequence
#    (matched per-&-atom) + panel gene + rare. Consequence is split on '&' so
#    unanticipated compound terms are not silently dropped. A panel-gene variant
#    that is ClinVar P/LP with >=1 review star is EXEMPT from the rarity and
#    consequence gates (MANE + panel still required), so founder alleles above the
#    AF ceiling and pathogenic non-coding/synonymous variants are not discarded
#    before any evidence arm runs.                                     [#1]
#  * Frequency threshold is MOI-aware: recessive genes tolerate a higher
#    gnomAD AF than dominant genes.                                    [#6]
#  * Inclusion (rescue) gate (OR): CADD>=25.3, AlphaMissense>=0.792, EVE path,
#    REVEL>=0.644, Pangolin>=0.5, ClinVar P/LP, PS1/PM5 (ClinVar AA match), LoF
#    (LOFTEE HC or high-impact truncating consequence), AR_hom (see below). Each
#    surviving row records which arm(s) fired in a `kept_by` column, using the
#    tokens CADD/AM/EVE/REVEL/Pangolin/ClinVar/PS1/PM5/LoF/AR_hom.   [#4,#8]
#  * Genotype-aware: proband zygosity / DP / GQ / allele-balance columns from
#    FORMAT; parent inheritance uses parental GT (carrier = non-ref), not mere
#    site presence.                                                [#5,#10]
#  * One MANE row per variant: a variant annotating onto >1 MANE transcript
#    (MANE Select + MANE Plus Clinical, or overlapping gene models) is collapsed to
#    a single row — prefer panel-primary, then MANE Select, then most evidence arms.
#    (--lookup consults still report every annotation.)                     [#3]
#  * Per-gene recessive logic (recessive-capable genes only): homozygous and
#    compound-het (trans where phaseable) flags. Purely dominant genes get no
#    HOM/CompHet flag, so two independent hets are not mislabeled comp-het. [#6,#7]
#  * Dual-inheritance genes (MOI lists BOTH AD and AR, e.g. "AD, AR") are treated
#    as DOMINANT for the carrier drop: a solitary het passes through as a candidate
#    (no recessive flag), while a genuine hom/comp-het still gets the recessive
#    flag. Pure AR/XLR genes use the recessive path below.                    [#5]
#  * Recessive carriers (DEFAULT = drop): a solitary het in a PURE recessive (AR/XLR)
#    gene that is not biallelic is DROPPED (a single het can't explain recessive disease;
#    carrier states are clinical noise). True comp-hets are unaffected — a gene with >=2
#    gate-passing hets is a biallelic CompHet and kept. OPT-IN --keep-ar-carriers /
#    KEEP_AR_CARRIERS=1 surfaces the STRONG such carriers (carrier-only tier: ClinVar
#    P/LP >=1*, HC-LoF, or >=2 strong predictors AM>=0.906/CADD>=28.1/EVE-path/REVEL>=
#    0.773, not Benign/LB; flagged flags=carrier-only) for second-hit hunts. [#1,#2,#6]
#  * ClinVar (fresh, via --custom), gnomAD nhomalt + FILTER surfaced as
#    columns.                                                      [#4,#7]
#  * ACMG SF secondary findings: the 81 ACMG SF v3.2 genes are ALWAYS scanned
#    (independent of the candidate panel) with a STRICTER gate — ClinVar P/LP
#    (>=1 star, frequency-agnostic) OR novel LOFTEE-HC OR >=2 strong computational
#    predictors (AM>=0.906, CADD>=28.1, EVE path, REVEL>=0.773); AR genes report
#    biallelic only. These appear in the SAME candidatos output flagged with
#    GDV=Incidental (Association/MOI from the ACMG table; kept_by = evidence tier).
#  * Automated ACMG/AMP classification (TRIAGE ONLY): per-row acmg_class +
#    acmg_criteria, combined per the categorical ACMG 2015 rules. PP3/BP4 use a
#    single CALIBRATED tool — AlphaMissense primary (Bergquist 2025), REVEL
#    fallback (Pejaver 2022) — graded Supporting/Moderate/Strong with a REVEL
#    direction-conflict veto, mapped to the 2015 strength tiers (BP4_Moderate ->
#    supporting-benign, as 2015 has no benign-Moderate). PP2: missense in a gene with
#    low benign-missense variation, from gnomAD v4.1.1 missense constraint (mis.oe <
#    0.6 on MANE, outliers excluded); counts independently of PP3, but suppressed when
#    BP4 fired (no gene-level pathogenic support for a benign-predicted variant). Other
#    criteria: PVS1, PS2/PM6, PM2/PM4, PP5, BA1/BS1/BS2/BP6/BP7.           [#2]
#  * Cohort recurrent-artifact filter (internal panel-of-normals): for a real
#    cohort auto-analyzed together (>= $COHORT_MIN probands), a candidate carried by
#    >= $COHORT_MAX_FRAC of samples AND absent from gnomAD (joint AC == 0) is
#    a systematic technical artifact (paralog/low-complexity mismapping) and is
#    dropped — both conditions required, so population bottlenecks / founder alleles
#    (which carry a gnomAD footprint) are preserved. OFF for single-variant, forced/
#    single-proband, and small runs. --keep-cohort-artifacts keeps+tags instead. [#11]
#  * QC / artifact flags (in the consolidated `flags` column, after the recessive
#    verdict): lowDP, lowGQ, AB_het/AB_hom, homopolymer, clinvar_conflict,
#    (indels, via samtools+reference), inh_lowqual, DN_unconfirmed.   [#6,#7]
#    NOTE: parent VCFs are variant-only, so de-novo cannot be confirmed from
#    parental reference depth — DN is flagged DN_unconfirmed by design.      [#6]
#  * Run summary printed per proband.                                  [#9]
#
# Splicing (Pangolin) two-pass bridge
# -----------------------------------
#  If <proband>.<panel>.pangolin.tsv is ABSENT, the script writes
#  <proband>.<panel>.pangolin_input.csv (the structural-pass candidates) and
#  stops for that proband. Run Pangolin (run_filtering.sh) to create the .tsv,
#  then re-run to produce <proband>.<panel>.candidatos.  (<panel> = panel
#  basename, so different gene lists produce side-by-side outputs.)
#
# Single-variant lookup mode
# --------------------------
#  Report EVERYTHING for one (or a few) variants with the panel / rarity /
#  consequence / rescue gates and the recessive-carrier drop all BYPASSED
#  (off-panel genes get Association/MOI/GDV = NA; kept_by lists whichever evidence
#  arms fire, else "none"). MANE-only unless --all-transcripts. Genotype columns
#  are blank (sites-only) and inheritance = NA.
#  Output is the TRANSPOSED, human-readable view only — one "field <TAB> value"
#  line per column (not the TSV candidatos table) — written to
#  Lookup.<coords>.<panel>.candidatos and echoed to stdout.
#  Two ways in:
#    -v/--variant '<v>'   resolve + annotate variant(s) from scratch (repeatable).
#    --lookup <vcf.gz>    analyze a pre-annotated *.germline.vep.vcf.gz directly.
#  <v> is dashed/colon GRCh38 coords (chr2-166199981-A-G | 2:166199981:A:G, also
#  the 5-field chr-start-end-ref-alt form), resolved offline, OR transcript-
#  qualified HGVS (ENST…:c.…), recoded to coords via the Ensembl REST API (only the
#  variant string is sent, never patient data).
#  -v builds a sites-only VCF and runs vep_annotate.sh, removing the annotated
#  VCF afterward unless --keep-vcf. The panel (Association/MOI/GDV columns) is the
#  default g4e-2026 unless overridden with -l/--list <genes>.
#############################################################################

# ── Tunable thresholds (single source of truth) ──
my $INPUT_GLOB = '*.germline.vep.vcf.gz';
my $FREQ_AD    = 0.01;   # max gnomAD AF (%) for dominant genes  (1 in 10,000)
my $FREQ_AR    = 1.0;    # max gnomAD AF (%) for recessive genes (carrier freq)
my $CADD_MIN   = 25.3;   # CADD PHRED rescue threshold
my $REVEL_MIN  = 0.644;  # REVEL rescue threshold (ClinGen PP3)
my $AM_MIN     = 0.792;  # AlphaMissense pathogenicity rescue threshold
my $SPLICE_MIN = 0.5;    # Pangolin |delta| splice rescue threshold

# ── Splice DISCOVERY probes [#5] ──
# typevar.txt has no bare `intron_variant` and no bare `synonymous_variant`, so a
# deep-intronic or exonic-synonymous splice-disrupting variant was dropped at Stage 1,
# never reached Pangolin, and the splice rescue arm could only ever UPGRADE a variant that
# was already whitelisted — never DISCOVER one. The classic pathogenic deep-intronic
# alleles (CFTR c.3718-2477C>T, USH2A c.7595-2144A>G) were structurally unreachable.
#
# Probes are scored by Pangolin and kept ONLY if the splice arm fires (>= $SPLICE_MIN).
# They are not candidates in their own right, so a probe that Pangolin scores low simply
# disappears — the probe set widens what can be FOUND without widening the table.
#
# VOLUME IS THE BINDING CONSTRAINT. Pangolin scores ~8 variants/s on one GPU, and a WGS
# proband carries ~155k rare intronic variants in panel genes. Two bounds keep it finite:
#   (a) distance from the exon boundary, read from the HGVSc offset; and
#   (b) a strict rarity ceiling applied to ALL genes — a pathogenic splice variant is rare
#       regardless of the gene's mode of inheritance, so the permissive $FREQ_AR carrier
#       ceiling (1%, meant for recessive coding candidates) is deliberately NOT used here.
# Set $SPLICE_PROBE = 0, or pass --no-splice-discovery, to restore the prior behaviour.
my $SPLICE_PROBE    = $ENV{NO_SPLICE_DISCOVERY} ? 0 : 1;
my $INTRON_MAX_DIST = 300;    # max |HGVSc intron offset| to probe (bp from exon boundary)
my $PROBE_FREQ_MAX  = 0.01;   # max gnomAD AF (%) for a probe, every gene

# Require the probe position to EXIST in the gnomAD resource (AN > 0). Default ON.
#
# MEASURED COST OF THIS SETTING. The custom VCF is gnomAD.joint.v4.1.**mane**, so intronic
# coverage is thin: on a real WGS proband the 300 bp window contains ~22,700 rare intronic
# variants, but only ~8 of them are gnomAD-covered. Requiring coverage therefore keeps the
# probe set almost free (8 probes/proband) — and keeps discovery confined to roughly the
# MANE footprint, so a classic deep-intronic allele (CFTR c.3718-2477C>T) is still out of
# reach. Turning it OFF (--probe-uncovered) probes all ~22,700: at ~8 variants/s that is
# ~48 min per proband of GPU time, and every rescued variant sits in territory where the
# rarity gate cannot work. PM2 is withheld on such rows (see acmg_classify) so the run at
# least does not invent a pathogenic criterion for them.
#
# The real unlock is rebuilding the custom gnomAD VCF with genome-wide coverage; then this
# flag stops mattering and the rarity gate works everywhere.
my $PROBE_REQUIRE_GNOMAD = $ENV{PROBE_UNCOVERED} ? 0 : 1;

# ── Cohort recurrent-artifact filter (internal panel-of-normals) [#11] ──
# Systematic technical artifacts — reference/mapping errors in paralog-rich or
# low-complexity genes (e.g. SYNE1, KMT2C) — recur across a large fraction of a
# cohort yet are ABSENT from gnomAD. That combination is the discriminator: a
# population bottleneck or an under-represented ancestry CANNOT produce it, because
# a real founder allele frequent enough to reach a quarter of the cohort would
# leave a footprint in gnomAD's large Admixed-American sample. So a variant is
# dropped ONLY when it is BOTH cohort-recurrent (>= $COHORT_MAX_FRAC of samples) AND
# wholly absent from gnomAD (joint AC == 0). Neither condition alone drops anything.
# Absence is tested as AC == 0 rather than against a frequency ceiling on purpose: a
# ceiling at or above $FREQ_AD is already implied by the Stage-1 rarity gate for
# dominant genes, which would leave recurrence as the only effective condition and
# strip the founder-safety the second condition exists to provide. Requiring a
# literal zero keeps that protection meaningful for every mode of inheritance —
# any gnomAD footprint at all, however small, spares the variant.
# Activates ONLY for a real cohort auto-analyzed together (>= $COHORT_MIN probands,
# no forced selection): single-variant (-v/--lookup) and single/forced-proband
# (--proband) runs never activate it, and trios/duos/small runs are untouched
# (a per-variant or per-proband consult has no cohort to compare against). Logged;
# --keep-cohort-artifacts (env KEEP_COHORT_ARTIFACTS=1) keeps them instead, tagged
# flags=cohort_artifact — for a founder-enriched cohort, review the drop log, as a
# genuinely private founder allele would surface there.
# THRESHOLDS (revised 2026-08). The filter previously needed >= 10 probands, which no
# real internal batch reaches — batches run 6-9 samples — so it had never once fired on
# a clinical run, and the recurrent KMT2C/SYNE1 mismapping artifacts reached every
# delivered table. Two changes make it work at the batch sizes actually used:
#
#   1. $COHORT_MIN 10 -> 5 samples. Five unrelated genomes are enough to tell a
#      systematic artifact from a private allele when gnomAD-absence is also required.
#   2. A new ABSOLUTE floor, $COHORT_MIN_CARRIERS, on top of the fraction. At N=8 the
#      25% fraction alone means 2 carriers, and two unrelated Chilean probands sharing
#      a gnomAD-absent allele is an ordinary founder/relatedness event, not evidence of
#      a technical artifact. Requiring >= 3 carriers AND >= 25% keeps both ends honest:
#      the floor binds on small batches, the fraction binds on large cohorts (at N=56,
#      25% = 14 carriers, well above the floor).
#
# Measured on batch4 (8 singleton probands), the carrier distribution is bimodal with a
# clean gap: the four artifacts sit at 5,6,7,8 of 8 carriers (63-100%) and the next most
# recurrent candidate is 2 of 8 (25%). Both new thresholds fall inside that gap.
my $COHORT_MIN          = 5;     # min samples to activate the filter (below -> OFF)
my $COHORT_MAX_FRAC     = 0.25;  # carried by >= this fraction of the cohort -> "recurrent"
my $COHORT_MIN_CARRIERS = 3;     # ...AND by at least this many samples, whatever the fraction

# Ensembl REST endpoint for HGVS->coordinate recoding in -v variant mode
# (override with $ENSEMBL_REST; point at a private mirror on an air-gapped host).
my $ENSEMBL_REST = $ENV{ENSEMBL_REST} // 'https://rest.ensembl.org';

# High-impact loss-of-function consequences (LoF rescue arm; see inclusion gate).
my %LOF_CONS = map { $_ => 1 }
    qw(frameshift_variant stop_gained splice_acceptor_variant splice_donor_variant start_lost);

# ── ACMG SF secondary findings (always evaluated; stricter than candidates) ──
# Emitted into the SAME candidatos output, flagged GDV=Incidental.
my $ACMG_FILE   = 'acmg_sf_v3.2.txt';
my $SF_FREQ_MAX = 0.5;     # max gnomAD AF (%) for the NOVEL SF tiers (LoF/computational)
my $SF_AM       = 0.906;   # AlphaMissense pathogenicity (strong)
my $SF_CADD     = 28.1;    # CADD PHRED (strong)
my $SF_REVEL    = 0.773;   # REVEL (ClinGen PP3_moderate)

# ── QC / artifact flags [#7] and parental-quality de-novo confidence [#6] ──
my $QC_MIN_DP   = 15;      # depth below this -> lowDP
my $QC_MIN_GQ   = 20;      # genotype quality below this -> lowGQ
my $REF_FASTA   = $ENV{REF_FASTA} // $ENV{PANGOLIN_FASTA} // '';  # set by site.sh; blank => check skipped
my $HAVE_REF    = -e "$REF_FASTA.fai";   # samtools-indexed reference for homopolymer check

# ── Automated ACMG/AMP classification (InterVar-style, triage only) [#2] ──
my $PM2_AC_MAX  = 1;       # gnomAD AC at/below -> PM2 (absent=0 or singleton=1)
# PM2 evidence strength. ACMG 2015 lists PM2 as Moderate; ClinGen SVI (2020) recommends
# downgrading it to SUPPORTING for rare disease, because absence from a population database
# is weak evidence on its own and is the single most over-applied criterion.
#
# DELIBERATELY LEFT AT 'moderate'. The SVI downgrade is only coherent inside the framework
# it was calibrated in — the ClinGen/Tavtigian Bayesian POINTS system, where PVS1=8 and
# PM2_Supporting=1 sum to 9 points and still reach Likely pathogenic. This pipeline combines
# CATEGORICALLY (ACMG 2015 Table 5), which has no "PVS1 + 1 supporting" pathway at all, so
# the downgrade there silently demotes every gnomAD-absent nonsense/frameshift variant in a
# disease gene to VUS. Measured on batch4: KCNT1, CUX2, RELN and HCN2 LoF calls all dropped
# to VUS. That is an artefact of mixing two frameworks, not a more conservative reading.
#
# Setting this to 'supporting' is therefore only correct once the combining step is replaced
# with a points-based one. Until then 'moderate' is the self-consistent choice.
# See the README limitations section.
my $PM2_STRENGTH = 'moderate';     # 'moderate' (ACMG 2015) | 'supporting' (ClinGen SVI 2020)
my $BS1_FREQ    = 1.0;     # gnomAD AF (%) at/above -> BS1 (too common for rare disease)
my $BA1_FREQ    = 5.0;     # gnomAD AF (%) at/above -> BA1 (benign standalone)
my $BS2_NHOM    = 10;      # gnomAD homozygotes at/above -> BS2
my $BP4_REVEL   = 0.290;   # REVEL at/below -> BP4 (computational benign, ClinGen)
my $PP2_MIS_OE  = 0.6;     # gnomAD v4.1.1 missense o/e below this -> gene is missense-
                           # constrained; a missense variant there earns PP2 (see acmg_classify)
my $CONSTRAINT_FILE = 'gnomad-mis-constraint.txt';   # gene -> mis.oe / mis.z / flags

#############################################################################
# Reference hashes
#############################################################################

# Cohort recurrent-artifact self-test (synthetic VCFs; no annotation stack needed).
#   perl filtering_r.pl --selftest-cohort
run_cohort_selftest() if grep { $_ eq '--selftest-cohort' } @ARGV;

# ── Where the tracked reference data lives ──
# Resolve each reference file from the CURRENT directory first (so a run directory can
# drop in its own panel), then from the script's own directory. Without this every
# reference file was opened by bare relative name, so the documented $WORKDIR override
# only worked if you symlinked the whole repo into the run directory — which is exactly
# what the per-sample batch directories had been doing.
my $CF_REPO = $0;
$CF_REPO = ($CF_REPO =~ s{/[^/]+$}{}r);
$CF_REPO = "." if $CF_REPO eq $0 || $CF_REPO eq "";
sub ref_file {
    my ($f) = @_;
    return $f if -e $f;
    return "$CF_REPO/$f" if -e "$CF_REPO/$f";
    return $f;                      # let the caller's open() report the error
}

open MANE, "<", ref_file("mane-plus-clinical-names.txt") or die "mane: $!";
my %mane;
while (my $m = <MANE>) { chomp $m; my ($a,$b) = split /\t/, $m; $mane{$a} = $b; }
close MANE;
print "hash mane, listo!\n";

# gnomAD v4.1.1 missense constraint (MANE Select): gene -> mis.oe, for the ACMG PP2
# criterion. Genes flagged as missense-constraint outliers (outlier_mis / no_exp_mis)
# are skipped so they can never earn PP2. Optional file — absent -> PP2 simply never fires.
# gnomAD v4.1.1 predates several HGNC symbol changes, so a current panel/VEP symbol
# can miss its own constraint record and silently forfeit PP2. These five are confirmed
# renames where the OLD symbol is present in the constraint file (verified by lookup);
# every other panel gene absent from the file is genuinely absent (mitochondrial, snRNA/
# snoRNA, or not MANE Select in v4.1.1) and has no missense constraint to inherit.
my %GENE_ALIAS = (
    GBA1  => 'GBA',        # HGNC 2022
    BMAL1 => 'ARNTL',      # HGNC 2023
    AFG2A => 'SPATA5',     # HGNC 2022
    AFG2B => 'SPATA5L1',   # HGNC 2022
    BLTP1 => 'KIAA1109',   # HGNC 2023
);

my %mis_oe;
if (open my $mc, "<", ref_file($CONSTRAINT_FILE)) {
    while (my $l = <$mc>) {
        next if $l =~ /^#/ || $l !~ /\S/;
        chomp $l;
        my ($g,$oe,$z,$flags) = split /\t/, $l;
        next unless defined $oe && $oe ne "";
        next if defined $flags && $flags =~ /outlier_mis|no_exp_mis/;
        $mis_oe{$g} = $oe + 0;
    }
    close $mc;
    printf "missense constraint (PP2): %d MANE genes loaded (mis.oe < %s -> constrained)\n",
           scalar keys %mis_oe, $PP2_MIS_OE;
} else {
    warn "NOTE: $CONSTRAINT_FILE not found — ACMG PP2 disabled (missense constraint unavailable)\n";
}

# Missense o/e for a gene symbol, falling back to its pre-rename symbol. Returns
# undef when the gene has no constraint record (PP2 then simply never fires).
sub mis_oe_for {
    my ($g) = @_;
    return undef unless defined $g && $g ne "";
    return $mis_oe{$g} if exists $mis_oe{$g};
    my $alias = $GENE_ALIAS{$g};
    return (defined $alias && exists $mis_oe{$alias}) ? $mis_oe{$alias} : undef;
}

# ── Argument parsing ──
#   --proband NAME   (repeatable) force NAME as a proband, overriding filename-
#                    based auto-discovery. NAME must have <NAME>.germline.vep.vcf.gz.
#   -l/--list FILE   candidate-gene panel override (replaces default g4e-2026).
#                    The ONLY way to set the panel (no positional form). Applies
#                    to normal runs and to -v/--lookup variant consults alike.
#
# Gene panel: gene -> "Association \t MOI \t GDV".
#   Default source = g4e-2026.txt (4 columns). A genes-of-interest file given via
#   -l/--list (one gene symbol per line; plain symbols -> Association/MOI/GDV =
#   "NA", or full 4-column g4e format) overrides it. '#' comments/blanks skipped.
my (@force_probands, $GENES_FILE, $LOOKUP_FILE, @VARIANTS);
my ($LOOKUP, $ALL_TX, $KEEP_VCF, $NO_SPLICE) = (0, 0, 0, 0);
# Carrier opt-in. By DEFAULT a solitary het carrier of a pure recessive (AR/XLR) gene is
# DROPPED (biallelic-only; carrier states are clinical noise, and true comp-hets are kept
# independently via the CompHet flag). Set KEEP_AR_CARRIERS=1 or --keep-ar-carriers to
# SURFACE the strong such carriers (carrier-only tier: strong-evidence & not benign,
# flagged flags=carrier-only; see the [#1,#2] block) — e.g. to chase a possible
# missed second hit (deep-intronic, CNV) in a targeted investigation.
my $KEEP_AR_CARRIERS = $ENV{KEEP_AR_CARRIERS} ? 1 : 0;
# Keep (don't drop) cohort recurrent-artifact variants, tagging them flags=
# cohort_artifact instead. Env KEEP_COHORT_ARTIFACTS=1 or CLI --keep-cohort-artifacts.
my $KEEP_COHORT_ARTIFACTS = $ENV{KEEP_COHORT_ARTIFACTS} ? 1 : 0;
while (@ARGV) {
    my $a = shift @ARGV;
    if    ($a eq '--keep-ar-carriers')      { $KEEP_AR_CARRIERS = 1; }
    elsif ($a eq '--keep-cohort-artifacts') { $KEEP_COHORT_ARTIFACTS = 1; }
    elsif ($a eq '--proband' || $a eq '-p') { push @force_probands, (shift @ARGV // ''); }
    elsif ($a =~ /^--proband=(.+)$/)        { push @force_probands, $1; }
    elsif ($a eq '--lookup')                { $LOOKUP = 1; $LOOKUP_FILE = (shift @ARGV // ''); }
    elsif ($a =~ /^--lookup=(.+)$/)         { $LOOKUP = 1; $LOOKUP_FILE = $1; }
    elsif ($a eq '--variant' || $a eq '-v') { push @VARIANTS, (shift @ARGV // ''); }
    elsif ($a =~ /^--variant=(.+)$/)        { push @VARIANTS, $1; }
    elsif ($a eq '--list'  || $a eq '-l')   { $GENES_FILE = (shift @ARGV // ''); }   # candidate-gene panel override
    elsif ($a =~ /^--list=(.+)$/)           { $GENES_FILE = $1; }
    elsif ($a eq '--all-transcripts')       { $ALL_TX = 1; }
    elsif ($a eq '--keep-vcf')              { $KEEP_VCF = 1; }
    elsif ($a eq '--no-splice')             { $NO_SPLICE = 1; }   # -v: skip Pangolin splice scoring
    elsif ($a eq '--no-splice-discovery')   { $SPLICE_PROBE = 0; } # skip the intronic/synonymous probe set
    elsif ($a eq '--probe-uncovered')       { $PROBE_REQUIRE_GNOMAD = 0; }  # probe outside the gnomAD footprint (SLOW)
    else {
        die "unknown argument: '$a'\n".
            "  set the candidate-gene panel with  -l/--list FILE  (no positional form);\n".
            "  a single variant with  -v/--variant '<v>' , a pre-annotated VCF with  --lookup FILE;\n".
            "  keep recessive carriers with  --keep-ar-carriers , cohort artifacts with  --keep-cohort-artifacts;\n".
            "  skip the splice-discovery probe set with  --no-splice-discovery ,\n".
            "  or widen it past the gnomAD footprint (SLOW) with  --probe-uncovered .\n";
    }
}
my $PANEL = (defined $GENES_FILE && $GENES_FILE ne "") ? $GENES_FILE : "g4e-2026.txt";
my $custom_panel = (defined $GENES_FILE && $GENES_FILE ne "") ? 1 : 0;

# Output tag = panel basename without extension, minus any trailing year suffix
# (e.g. g4e-2026.txt -> g4e, Hyperparathyroidism.txt -> Hyperparathyroidism).
# All per-run outputs are namespaced by it so different panels don't overwrite.
my $PANEL_TAG = $PANEL;
$PANEL_TAG =~ s{.*/}{};
$PANEL_TAG =~ s/\.[^.]+$//;
$PANEL_TAG =~ s/-\d{4}$//;               # drop trailing year suffix (g4e-2026 -> g4e)

open PANEL, "<", ref_file($PANEL) or die "gene panel '$PANEL': $!";
my %epigenes;
while (my $g = <PANEL>) {
    $g =~ s/\r?\n$//;
    $g =~ s/^\s+|\s+$//g;
    next if $g eq "" || $g =~ /^#/;
    my @f = split /\t/, $g;
    my $sym = $f[0];
    $epigenes{$sym} = (@f >= 4) ? join("\t", $f[1], $f[2], $f[3]) : "NA\tNA\tNA";
}
close PANEL;
printf "gene panel: %s (%d genes%s)\n", $PANEL, scalar(keys %epigenes),
       $custom_panel ? ", custom — missing Association/MOI/GDV = NA" : "";

# PP2 coverage. A panel gene with no constraint record can never earn PP2, and the
# failure is otherwise invisible — report it once at startup instead of letting the
# criterion go quietly missing for part of the panel.
if (%mis_oe) {
    my @no_constraint = sort grep { !defined mis_oe_for($_) } keys %epigenes;
    printf "PP2 coverage: %d/%d panel genes have gnomAD missense constraint%s\n",
           scalar(keys %epigenes) - scalar(@no_constraint), scalar(keys %epigenes),
           @no_constraint ? " — no record for: ".join(" ", @no_constraint) : "";
}

# Consequence whitelist (atomic terms recommended; compound entries harmless).
open VAR, "<", ref_file("typevar.txt") or die "typevar: $!";
my %varfilter;
while (my $t = <VAR>) { chomp $t; my ($c,$d) = split /\t/, $t; $varfilter{$c} = $d//""; }
close VAR;
print "hash var, listo!\n";

# ACMG SF v3.2 genes: gene -> "condition \t MOI \t report_category". Always loaded
# (independent of the candidate panel). Non-fatal if absent.
my %acmg;
if (open my $afh, "<", ref_file($ACMG_FILE)) {
    while (my $g = <$afh>) {
        chomp $g; next if $g =~ /^#/ || $g !~ /\S/;
        my ($sym,$cond,$moi,$cat) = split /\t/, $g;
        $acmg{$sym} = join("\t", $cond//"", $moi//"AD", $cat//"ALL_PLP");
    }
    close $afh;
    printf "ACMG SF: %d genes (secondary findings, always evaluated)\n", scalar keys %acmg;
} else {
    warn "WARN: $ACMG_FILE not found — secondary findings (Incidental) disabled\n";
}

# ── ClinVar amino-acid evidence for PS1/PM5 (optional; graceful if absent) ──
# Built from the MANE-missense split; override the directory with $CLINVAR_AA_DIR.
my $CLINVAR_AA_DIR = $ENV{CLINVAR_AA_DIR} // '';   # set CLINVAR_AA_DIR in site.env (README section 0.6)
my $PLP_resid = load_clinvar_aa("$CLINVAR_AA_DIR/clinvar.MANE_missense.PLP.tsv");
my $BLB_resid = load_clinvar_aa("$CLINVAR_AA_DIR/clinvar.MANE_missense.BLB.tsv");
my $CLINVAR_AA_ON = (keys %$PLP_resid) ? 1 : 0;
if ($CLINVAR_AA_ON) {
    printf "ClinVar AA evidence (PS1/PM5): %d P/LP residues, %d B/LB residues\n",
        scalar(keys %$PLP_resid), scalar(keys %$BLB_resid);
} else {
    warn "WARN: ClinVar missense P/LP resource not found"
       . ($CLINVAR_AA_DIR ? " under $CLINVAR_AA_DIR" : " (CLINVAR_AA_DIR unset)")
       . " — PS1/PM5 disabled. See README section 0.6.\n";
}

#############################################################################
# Helpers
#############################################################################

sub open_vcf {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) { open($fh,"-|","gzip","-dc",$file) or die "gzip $file: $!"; }
    else                  { open($fh,"<",$file)              or die "$file: $!"; }
    return $fh;
}

# CSQ field name -> column index, from the ##INFO=<ID=CSQ ... Format: ...> header.
sub csq_columns {
    my ($file) = @_;
    my $fh = open_vcf($file);
    my %col;
    while (my $line = <$fh>) {
        last if $line =~ /^#CHROM/;
        next unless $line =~ /ID=CSQ/;
        if ($line =~ /Format:\s*([^"]+)"/) {
            my @n = split /\|/, $1;
            $col{$n[$_]} = $_ for 0 .. $#n;
        }
        last if %col;
    }
    close $fh;
    return \%col;
}

# Resolve a logical field to a CSQ index: exact name(s) first, then regex.
sub resolve {
    my ($col, @cand) = @_;
    for my $name (@cand) { return $col->{$name} if exists $col->{$name}; }
    for my $pat (@cand) {
        for my $name (keys %$col) { return $col->{$name} if $name =~ /$pat/i; }
    }
    return undef;
}

sub field {
    my ($row,$i) = @_;
    return "" unless defined $i;
    my $v = $row->[$i];
    return defined($v) ? $v : "";
}

# Parse one sample's FORMAT:SAMPLE -> (GT, DP, GQ, AD_ref, AD_alt).
sub parse_call {
    my ($fmt,$smp) = @_;
    return ("","","","","") unless defined $fmt && defined $smp;
    my @k = split /:/, $fmt;
    my @v = split /:/, $smp;
    my %h; @h{@k} = @v;
    my $gt = defined $h{GT} ? $h{GT} : "";
    my $dp = defined $h{DP} ? $h{DP} : "";
    my $gq = defined $h{GQ} ? $h{GQ} : "";
    my ($ar,$aa) = ("","");
    if (defined $h{AD} && $h{AD} ne "" && $h{AD} ne ".") {
        my @ad = split /,/, $h{AD};
        ($ar,$aa) = ($ad[0]//"", $ad[1]//"");
    }
    return ($gt,$dp,$gq,$ar,$aa);
}

# Zygosity from a GT string: hom (alt/alt), het (ref/alt), hem (haploid alt), ref,
# or "" (no-call).
# HAPLOID CALLS: DRAGEN (and GATK with -ploidy 1) emit a single-allele GT — "1" or
# "0" — for non-PAR chrX/chrY in males and for chrM. Requiring two alleles made those
# calls return "", which silently cost them their zygosity (no HOM flag, no AR_hom
# rescue, no AB_hom QC) and, worse, made a hemizygous parent invisible to load_parent
# so an INHERITED X-linked variant was reported as de novo. A hemizygous call is a
# complete genotype, not half a het: it is reported as "hem" and treated as biallelic-
# equivalent everywhere the recessive logic asks "is this genotype sufficient?".
sub zygosity {
    my ($gt) = @_;
    return "" unless defined $gt && $gt ne "";
    my @a = split /[\/|]/, $gt;
    return "" unless @a && !grep { $_ eq "." || $_ eq "" } @a;
    my $n1 = grep { $_ eq "1" } @a;
    return $n1 ? "hem" : "ref" if @a == 1;          # haploid (non-PAR X/Y male, chrM)
    return $n1 >= 2 ? "hom" : $n1 == 1 ? "het" : "ref";
}

# Is this genotype a complete (non-carrier) recessive genotype on its own?
# Homozygous alt, or hemizygous alt on a haploid contig — both leave no second
# wild-type allele, so neither is a "carrier" state.
sub zyg_biallelic { my $z = shift // ""; return ($z eq "hom" || $z eq "hem") ? 1 : 0; }

# Parent carrier map: chr-pos-ref-alt -> "gt:dp:gq" if the parent carries the ALT
# (GT contains a '1'); records with 0/0 or no-call are NOT carriers. The DP/GQ are
# kept so inherited calls can report parental call quality [#6].
sub load_parent {
    my ($file) = @_;
    my %carry;
    return \%carry unless defined $file && -e $file;
    my $fh = open_vcf($file);
    while (my $line = <$fh>) {
        next if $line =~ /^#/;
        chomp $line;
        my @c = split /\t/, $line;
        my ($chr,$pos,$ref,$alt,$fmt,$smp) = @c[0,1,3,4,8,9];
        next if $alt =~ /,/;                       # should be pre-split
        my ($gt,$dp,$gq) = parse_call($fmt,$smp);
        my $z = zygosity($gt);
        # "hem" counts: a hemizygous father IS a carrier of the allele he transmits.
        $carry{"$chr-$pos-$ref-$alt"} = "$gt:$dp:$gq" if $z eq "het" || $z eq "hom" || $z eq "hem";
    }
    close $fh;
    return \%carry;
}

# Pangolin score map: "chr-pos-ref-alt <tab> score".
sub load_scores {
    my ($file) = @_;
    my %s;
    open(my $fh,"<",$file) or die "$file: $!";
    while (my $l = <$fh>) { chomp $l; my ($id,$v) = split /\t/, $l; $s{$id} = $v if defined $v; }
    close $fh;
    return \%s;
}

# ── Cohort recurrent-artifact tally (internal panel-of-normals) [#11] ──
# Scan every cohort VCF once (GENOTYPES ONLY — no CSQ parse, so it is cheap) and
# count, per chr-pos-ref-alt, how many distinct samples carry the ALT plus their
# zygosity breakdown. Returns (\%carriers, \%hom, \%het, $n_samples). Sites-only
# files (no sample column) contribute no carriers. Used only for large cohorts
# (see $COHORT_MIN); never called in lookup / single-proband runs.
sub build_cohort_tally {
    my ($files) = @_;
    my (%ac, %hom, %het);
    my $n = 0;
    for my $f (@$files) {
        next unless defined $f && -e $f;
        $n++;
        my $fh = open_vcf($f);
        while (my $line = <$fh>) {
            next if $line =~ /^#/;
            chomp $line;                           # else the last field keeps its newline
            my @c = split /\t/, $line;
            next unless @c >= 10;                  # need FORMAT + >=1 sample column
            my ($chr,$pos,$ref,$alt,$fmt,$smp) = @c[0,1,3,4,8,9];
            next if $alt =~ /,/;                    # pre-split; skip residual multiallelic
            my ($gt) = parse_call($fmt,$smp);
            my $z = zygosity($gt);
            next unless $z eq "het" || $z eq "hom" || $z eq "hem";
            my $id = "$chr-$pos-$ref-$alt";
            $ac{$id}++;
            zyg_biallelic($z) ? $hom{$id}++ : $het{$id}++;
        }
        close $fh;
    }
    return (\%ac, \%hom, \%het, $n);
}

# Pure artifact decision: carried by >= $COHORT_MAX_FRAC of the cohort AND absent
# from gnomAD (joint AC == 0), with the cohort large enough (>= $COHORT_MIN) to be
# meaningful. Kept pure (no I/O) so run_cohort_selftest can exercise it directly.
sub cohort_artifact_call {
    my ($carriers, $n, $gnomad_ac) = @_;
    return 0 unless defined $n && $n >= $COHORT_MIN;
    $carriers  ||= 0;
    $gnomad_ac ||= 0;
    return 0 if $carriers < $COHORT_MIN_CARRIERS;         # absolute floor (small batches)
    return (($carriers / $n) >= $COHORT_MAX_FRAC && $gnomad_ac <= 0) ? 1 : 0;
}

# Distance into an intron from the nearest exon boundary, taken from the HGVSc offset
# (c.1234+56A>G -> 56; c.1235-30A>G -> 30). A range keeps its closest endpoint
# (c.100+5_100+12del -> 5). Returns undef when the annotation carries no offset, which is
# the case for anything that is not intronic — so callers must also check the consequence.
sub intron_offset {
    my ($hgvsc) = @_;
    return undef unless defined $hgvsc && $hgvsc =~ /c\./;
    (my $c = $hgvsc) =~ s/^.*c\.//;
    my $min;
    while ($c =~ /[+-](\d+)/g) { $min = $1 if !defined $min || $1 < $min; }
    return $min;
}

# Is a consequence whitelisted? Pass if ANY '&'-separated atom is in the list.
sub consequence_ok {
    my ($csq) = @_;
    for my $atom (split /&/, $csq) { return 1 if exists $varfilter{$atom}; }
    return 0;
}

# ClinVar classification helpers (on a CLNSIG-style string).
sub clinvar_pathogenic {
    my ($s) = @_;
    return 0 unless defined $s && $s ne "";
    return 0 if $s =~ /conflicting/i;
    return 0 if $s =~ /benign/i;
    return ($s =~ /pathogenic/i) ? 1 : 0;
}
sub clinvar_benign {
    my ($s) = @_;
    return 0 unless defined $s && $s ne "";
    return ($s =~ /benign/i && $s !~ /pathogenic/i) ? 1 : 0;
}
# ClinVar review status -> star count (0-4). Handles both the VEP CLNREVSTAT
# underscore form ("criteria_provided,_single_submitter") and the ClinVar TSV
# space form ("criteria provided, single submitter"). The "no assertion criteria
# provided" string contains "criteria provided" — guarded explicitly to 0 stars.
sub clinvar_stars {
    my ($s) = @_;
    return 0 unless defined $s && $s ne "";
    my $t = lc $s; $t =~ tr/_/ /;
    return 4 if $t =~ /practice guideline/;
    return 3 if $t =~ /expert panel/;
    return 2 if $t =~ /multiple submitter/;
    return 0 if $t =~ /no assertion|no classification|no interpretation/;
    return 1 if $t =~ /single submitter|conflicting|criteria provided/;
    return 0;
}

# ── Mode-of-inheritance predicates over a panel/ACMG MOI string ──────────────
# A gene is "recessive-capable" if its MOI mentions AR/XLR/recessive, and
# "dominant-capable" if it mentions AD/XLD/dominant. Dual-inheritance genes
# (e.g. "AD, AR") satisfy BOTH — handled explicitly by the callers.
#
# PLAIN "XL" satisfies BOTH, deliberately. It is the Genes4Epilepsy vocabulary for
# an X-linked gene whose mechanism is not split into XLD/XLR — 72 of the 1078 g4e-2026
# genes, including CDKL5, MECP2, ARX, IQSEC2, PCDH19, DDX3X, ATRX, SLC6A8 and FLNA.
# Matching neither predicate (the previous behaviour) was an oversight, not a policy:
# those genes got no HOM/CompHet flag, never qualified for the AR_hom rescue, and were
# held to the strict dominant AF ceiling. Treating XL as dual-inheritance is the
# clinically safe reading — dominant-capable keeps a solitary het (these genes act
# dominantly in heterozygous females), recessive-capable earns a hemizygous/homozygous
# call its HEM/HOM flag and the AR_hom rescue. Same rule the "AD, AR" genes use.
# "XLR"/"XLD" keep their specific meaning: \bXL\b cannot match either.
sub moi_recessive { my $m = shift; return (defined $m && $m =~ /\bAR\b|XLR|\bXL\b|recessiv/i) ? 1 : 0; }
sub moi_dominant  { my $m = shift; return (defined $m && $m =~ /\bAD\b|XLD|\bXL\b|dominant/i)  ? 1 : 0; }

# ── Carrier-only tier gate ───────────────────────────────────────────────────
# Strong-evidence bar for SURFACING a solitary heterozygous carrier of a PURE
# recessive (AR/XLR) gene instead of dropping it (flagged flags=carrier-only).
# Mirrors the ACMG-SF strong gate: ClinVar P/LP (>=1 star) OR LOFTEE HC OR
# >=2 strong computational predictors (AM/CADD/EVE/REVEL). Takes a row data hashref.
sub carrier_strong_evidence {
    my ($d) = @_;
    return 1 if clinvar_pathogenic($d->{clinvar_sig}) && ($d->{clinvar_stars} // 0) >= 1;
    return 1 if ($d->{loftee} // "") eq "HC";
    my $n = 0;
    $n++ if defined $d->{am_score} && $d->{am_score} ne "" && $d->{am_score} >= $SF_AM;
    $n++ if defined $d->{cadd}     && $d->{cadd}     ne "" && $d->{cadd}     >= $SF_CADD;
    $n++ if ($d->{eve_class} // "") =~ /athogenic/;
    $n++ if defined $d->{revel}    && $d->{revel}    ne "" && $d->{revel}    >= $SF_REVEL;
    return $n >= 2 ? 1 : 0;
}

# A candidate whose only classification is Benign / Likely-benign (ACMG auto-class
# or a non-conflicting ClinVar B/LB that is not itself P/LP). Excluded from the
# carrier-only tier — a benign carrier is noise, not a candidate.
sub is_benign_class {
    my ($d) = @_;
    return 1 if lc($d->{acmg_class} // "") =~ /benign/;
    my $cs = lc($d->{clinvar_sig} // "");
    return 1 if $cs =~ /benign/ && $cs !~ /conflict/
             && !(clinvar_pathogenic($d->{clinvar_sig}) && ($d->{clinvar_stars} // 0) >= 1);
    return 0;
}

# Lexicographic ">" over two equal-length numeric preference arrays (used by the
# one-row-per-variant dedup to pick the row to keep).
sub _key_gt {
    my ($a, $b) = @_;
    for my $k (0 .. $#$a) {
        return 1 if $a->[$k] > $b->[$k];
        return 0 if $a->[$k] < $b->[$k];
    }
    return 0;
}

# ClinVar amino-acid-level evidence for PS1/PM5, built from the MANE-missense
# split (clinvar.MANE_missense.{PLP,BLB}.tsv). Returns a hashref keyed by
# "GENE\tAApos\tRefAA" -> { AltAA => { chr-pos-ref-alt => stars } }, so a candidate's
# residue can be checked for the SAME change (PS1) or a DIFFERENT pathogenic change
# (PM5). The innermost level records WHICH variant contributed each classification,
# which is what lets PS1 exclude the candidate's own ClinVar record (see aa_best_stars). A
# single-codon in-frame deletion of that residue also triggers PM5 (it is a
# different protein change at the same P/LP residue) — see the call site.
sub load_clinvar_aa {
    my ($file) = @_;
    my %resid;
    return \%resid unless defined $file && -e $file;
    open(my $fh, "<", $file) or do { warn "WARN: cannot read $file: $!\n"; return \%resid; };
    <$fh>;   # header
    while (my $l = <$fh>) {
        chomp $l;
        my @f = split /\t/, $l;
        # 1-based cols: 8 GeneSymbol, 11 ReviewStatus, 23 BB_AApos, 24 BB_RefAA, 28 AltAA
        my ($gene,$rev,$aapos,$refAA,$altAA) = @f[7,10,22,23,27];
        next unless defined $gene && defined $aapos && defined $refAA && defined $altAA;
        next if $gene eq "" || $aapos eq "" || $refAA eq "" || $altAA eq "" || $refAA eq $altAA;
        my $st = clinvar_stars($rev);
        my $k  = "$gene\t$aapos\t$refAA";
        # Keyed per SOURCE VARIANT (cols 3-6: Chr, PositionVCF, Ref, Alt), not just per
        # amino-acid change. PS1 requires a PREVIOUSLY established variant, so the record
        # belonging to the variant under classification has to be identifiable and
        # excluded — otherwise a variant that is itself ClinVar P/LP matches its own
        # submission and earns PS1 on top of PP5 from that one record.
        my ($vchr,$vpos,$vref,$valt) = @f[2,3,4,5];
        my $vid = (defined $vchr && defined $vpos && defined $vref && defined $valt)
                ? "$vchr-$vpos-$vref-$valt" : "";
        my $cur = $resid{$k}{$altAA}{$vid};
        $resid{$k}{$altAA}{$vid} = $st if !defined $cur || $st > $cur;
    }
    close $fh;
    return \%resid;
}

# Best review-star count among the ClinVar records carrying one amino-acid change,
# ignoring the record that IS the variant being classified ($self, "chr-pos-ref-alt";
# pass "" to consider every record). Returns 0 when nothing qualifies, so a variant
# whose only support is its own submission cannot earn PS1.
sub aa_best_stars {
    my ($by_variant, $self) = @_;
    return 0 unless ref $by_variant eq 'HASH';
    my $best = 0;
    for my $vid (keys %$by_variant) {
        next if defined $self && $self ne "" && $vid eq $self;
        my $st = $by_variant->{$vid};
        $best = $st if defined $st && $st > $best;
    }
    return $best;
}

# Is an INDEL in/adjacent to a homopolymer run (>=5)? Error-prone context. [#7]
my %hp_cache;
sub homopolymer_context {
    my ($chr,$pos,$ref,$alt) = @_;
    return 0 unless $HAVE_REF;
    return 0 if length($ref) == length($alt);           # SNV/MNV only flag indels
    my $key = "$chr-$pos";
    return $hp_cache{$key} if exists $hp_cache{$key};
    my $a = $pos - 12; $a = 1 if $a < 1;
    my $b = $pos + 12;
    my $seq = qx(samtools faidx "$REF_FASTA" "$chr:$a-$b" 2>/dev/null);
    $seq =~ s/^>.*\n//; $seq =~ s/\s+//g;
    my $hp = ($seq ne "" && $seq =~ /(.)\1{4,}/) ? 1 : 0;
    return $hp_cache{$key} = $hp;
}

# Calibrated PP3/BP4 thresholds — AlphaMissense (Bergquist et al., GIM 2025) and
# REVEL (Pejaver et al., AJHG 2022).
my %AMP = (
    am_pp3_strong=>0.990, am_pp3_mod=>0.906, am_pp3_supp=>0.792,   # AM PP3 (no BP4 strong)
    am_bp4_mod  =>0.099, am_bp4_supp=>0.169,                       # AM BP4
    rv_pp3_strong=>0.932, rv_pp3_mod=>0.773, rv_pp3_supp=>0.644,   # REVEL PP3
    rv_bp4_strong=>0.016, rv_bp4_mod=>0.183, rv_bp4_supp=>0.290,   # REVEL BP4
);

# Automated ACMG/AMP classification (TRIAGE ONLY — not a final clinical call).
# Criteria are combined per the categorical ACMG 2015 rules. PP3/BP4 come from a
# single CALIBRATED tool — AlphaMissense primary, REVEL fallback — graded
# Supporting/Moderate/Strong with a REVEL direction-conflict veto, then mapped to
# the 2015 strength tiers (PP3_Strong->strong, _Moderate->moderate, _Supporting->
# supporting; BP4_Strong->strong-benign, BP4_Moderate/_Supporting->supporting-
# benign, since the 2015 framework has no benign-Moderate tier). [#2]
sub acmg_classify {
    my (%v) = @_;
    my (@P,@B);

    # Pathogenic criteria
    my $pvs1 = ($v{loftee} eq "HC" || ($v{lof_type} && $v{loftee} ne "LC")) ? 1 : 0;
    push @P, "PVS1" if $pvs1;
    if ($v{inh} eq "DN") {                        # trio de novo (relatedness assumed)
        push @P, ($v{gt_clean} ? "PS2" : "PM6");
    } elsif ($v{inh} =~ m{^DN/} && $v{de_novo_mech}) { push @P, "PM6"; }
    # PM4 is evidence for a protein-length change; PVS1 already covers the loss-of-
    # function reading of the same event. VEP compound terms make them collide
    # (start_lost&inframe_deletion, frameshift_variant&stop_lost) because $lof_type is
    # matched per '&'-atom while this regex matches the whole string — two ACMG lines
    # from one protein-terminus effect, which pushes an LP call to Pathogenic.
    push @P, "PM4" if !$pvs1 && $v{consequence} =~ /inframe_(insertion|deletion)|stop_lost/;
    # PM2 at the configured strength (see $PM2_STRENGTH). Written as PM2_Supporting when
    # downgraded so the criteria string says which reading produced the class.
    # "Absent from gnomAD" is only assertable where gnomAD actually looked. A splice
    # probe rescued from outside the MANE-restricted resource has AN=0 because the position
    # is not IN the resource, not because the allele is unobserved — awarding PM2 there
    # would manufacture pathogenic evidence out of an annotation gap. Scoped to probe rows
    # so ordinary candidates, which sit inside the covered footprint, are unaffected.
    my $ac_assertable = !($v{probe} && ($v{an} // 0) <= 0);
    push @P, ($PM2_STRENGTH eq 'moderate' ? "PM2" : "PM2_Supporting")
        if $ac_assertable && $v{ac} ne "" && $v{ac} <= $PM2_AC_MAX;   # absent or singleton
    # PP5 requires >=1 review star, like every other ClinVar consumer in this file
    # (Stage-1 exemption, ACMG-SF tier, carrier tier, BP6). Without the gate a single
    # 0-star "no assertion criteria provided" submission — ~16% of the P/LP corpus —
    # supplied the criterion that lifts an LP call to Pathogenic.
    push @P, "PP5" if clinvar_pathogenic($v{clnsig}) && ($v{clnstar} // 0) >= 1;
    # PS1 (same AA change P/LP) or PM5 (different change, same residue P/LP), from
    # the ClinVar AA resource; conflicting matches are tagged but still counted (triage).
    push @P, $v{aa_crit} . ($v{aa_conflict} ? "(conflicting)" : "") if $v{aa_crit};

    # PP3 / BP4: single calibrated tool (AM primary, REVEL fallback), graded
    # Supporting/Moderate/Strong, with a REVEL direction-conflict veto.
    my ($am,$rv) = ($v{am_score}, $v{revel});
    my ($pp3,$bp4) = ("","");
    if ($am ne "") {                              # AlphaMissense primary
        $pp3 = ($am >= $AMP{am_pp3_strong}) ? "strong"
             : ($am >= $AMP{am_pp3_mod})    ? "moderate"
             : ($am >= $AMP{am_pp3_supp})   ? "supporting" : "";
        $bp4 = ($am <= $AMP{am_bp4_mod})    ? "moderate"
             : ($am <= $AMP{am_bp4_supp})   ? "supporting" : "";
        if ($rv ne "") {                          # REVEL direction-conflict veto
            $pp3 = "" if $pp3 && $rv <= $AMP{rv_bp4_supp};   # secondary calls benign
            $bp4 = "" if $bp4 && $rv >= $AMP{rv_pp3_supp};   # secondary calls pathogenic
        }
    } elsif ($rv ne "") {                          # REVEL fallback (AM absent)
        $pp3 = ($rv >= $AMP{rv_pp3_strong}) ? "strong"
             : ($rv >= $AMP{rv_pp3_mod})    ? "moderate"
             : ($rv >= $AMP{rv_pp3_supp})   ? "supporting" : "";
        $bp4 = ($rv <= $AMP{rv_bp4_strong}) ? "strong"
             : ($rv <= $AMP{rv_bp4_mod})    ? "moderate"
             : ($rv <= $AMP{rv_bp4_supp})   ? "supporting" : "";
    }
    push @P, "PP3_".ucfirst($pp3) if $pp3;
    push @B, "BP4_".ucfirst($bp4) if $bp4;

    # PM1: the variant falls in a PERv1 pathogenic-variant-enriched region naming
    # THIS gene (Perez-Palma et al., Genome Res 2020, whose stated application is PM1).
    # Two arms, both published:
    #   PERv1_direct  — enrichment computed on this gene. Graded by its own fold
    #                   enrichment at the Tavtigian-2018 calibration the paper cites:
    #                   >= 18.7 counts at Strong, else Moderate.
    #   PERv1_paralog — enrichment computed across the paralog family alignment and
    #                   assigned to every member, including members carrying none of
    #                   the underlying variants. This is the paper's headline arm
    #                   (1,252 genes vs 215) and the one its held-out de novo test
    #                   validated. Capped at Moderate upstream: the transfer step
    #                   costs a tier.
    # Restricted to missense and in-frame indels, which is what the region was
    # computed from -- a PER is a missense-burden statement, so it lends nothing to a
    # splice, synonymous or LoF call (and PVS1 already covers the last).
    # Deliberately NOT suppressed when BP4 fires: unlike PP2 (gene-level constraint),
    # PM1 here is regional evidence independently validated against held-out de novo
    # variants, so a benign computational prediction does not negate it.
    # $1 is captured BEFORE the consequence match: a second successful regex with no
    # capture groups clears it, which silently downgraded every PM1_Strong to PM1.
    my $pm1_st = ($v{pm1} // "") =~ /^(Strong|Moderate)$/ ? $1 : "";
    if ($pm1_st ne ""
        && $v{consequence} =~ /missense_variant|inframe_(?:insertion|deletion)/) {
        push @P, ($pm1_st eq "Strong" ? "PM1_Strong" : "PM1");
    }

    # PP2: missense in a gene with a low rate of benign missense variation, from
    # gnomAD v4.1.1 missense constraint (mis.oe < $PP2_MIS_OE on the MANE transcript;
    # constraint outliers already excluded at load). PP2 counts INDEPENDENTLY of PP3
    # (both are legitimate, separate ACMG lines — gene-level missense intolerance vs the
    # variant-level predictor — and ACMG 2015 permits combining them). It is still
    # suppressed when BP4 fired: a variant the calibrated tool predicts BENIGN must not
    # also collect gene-level pathogenic support (a genuine contradiction, not just
    # correlation). To let PP2 fire even alongside BP4, drop the "!$bp4" guard.
    push @P, "PP2" if $v{consequence} =~ /missense/
                   && defined $v{mis_oe} && $v{mis_oe} ne "" && $v{mis_oe} < $PP2_MIS_OE
                   && !$bp4;

    # Benign criteria
    push @B, "BA1" if $v{freq} >= $BA1_FREQ;
    push @B, "BS1" if $v{freq} >= $BS1_FREQ && $v{freq} < $BA1_FREQ;
    push @B, "BS2" if $v{nhom} ne "" && $v{nhom} >= $BS2_NHOM;
    push @B, "BP6" if clinvar_benign($v{clnsig}) && ($v{clnstar} // 0) >= 1;
    # BP7 requires POSITIVE evidence of no splice impact, so it fires only when a
    # Pangolin score actually exists. An unscored variant is unknown, not benign:
    # treating a missing score as 0 would assert benign-supporting evidence on
    # every synonymous variant in a run without Pangolin, biasing the class away
    # from the splice-active synonymous variants this pipeline exists to surface.
    push @B, "BP7" if $v{consequence} =~ /synonymous_variant/
                   && $v{pangolin} ne "" && $v{pangolin} < 0.2;

    # ── Combine per ACMG 2015. Count by strength tier; graded PP3/BP4 contribute
    #    at their tier (PP3_Strong->PS, _Moderate->PM, _Supporting->PP;
    #    BP4_Strong->BS, BP4_Moderate/_Supporting->BP). ──
    my $pvs = grep { $_ eq "PVS1" } @P;
    my $ps  = grep { /^PS\d/      } @P;            # PS1, PS2
    # PM2_Supporting must NOT be counted here: it starts with "PM2", so a bare /^PM\d/
    # would tally the downgraded criterion at Moderate anyway and silently undo it.
    # _Strong is excluded here for the same reason as _Supporting: "PM1_Strong"
    # starts with "PM1", so a bare /^PM\d/ would tally an upgraded criterion at
    # Moderate and silently undo the upgrade.
    my $pm  = grep { /^PM\d/ && !/_(?:Supporting|Strong)$/ } @P;   # PM2, PM4, PM5, PM6
    my $pp  = grep { $_ eq "PP2" || $_ eq "PP5" || /^PM\d_Supporting$/ } @P;
    $ps += grep { /^PM\d_Strong$/ } @P;            # PM1_Strong (PERv1 fold enrichment >= 18.7)
    $ps++ if $pp3 eq "strong";
    $pm++ if $pp3 eq "moderate";
    $pp++ if $pp3 eq "supporting";
    my $ba  = grep { $_ eq "BA1"  } @B;
    my $bs  = grep { /^BS\d/      } @B;            # BS1, BS2
    my $bp  = grep { /^BP[67]/    } @B;            # BP6, BP7
    $bs++ if $bp4 eq "strong";
    $bp++ if $bp4 eq "moderate" || $bp4 eq "supporting";

    my $path = ( ($pvs && ($ps >= 1 || $pm >= 2 || ($pm >= 1 && $pp >= 1) || $pp >= 2))
               || $ps >= 2
               || ($ps >= 1 && ($pm >= 3 || ($pm >= 2 && $pp >= 2) || ($pm >= 1 && $pp >= 4))) );
    my $lp   = ( ($pvs && $pm >= 1)
               || ($ps >= 1 && $pm >= 1)
               || ($ps >= 1 && $pp >= 2)
               || $pm >= 3 || ($pm >= 2 && $pp >= 2) || ($pm >= 1 && $pp >= 4) );
    my $ben  = ($ba || $bs >= 2);
    my $lb   = (($bs >= 1 && $bp >= 1) || $bp >= 2);

    my $pathy = $path || $lp;
    my $beny  = $ben  || $lb;
    my $class = ($pathy && $beny)  ? "Conflicting"
              :  $path             ? "Pathogenic"
              :  $lp               ? "Likely_pathogenic"
              :  $ben              ? "Benign"
              :  $lb               ? "Likely_benign"
              :                      "VUS";
    return ($class, join(",", @P, @B));
}

#############################################################################
# Single-variant entry point: resolve -> sites-only VCF -> vep_annotate.sh
#
# Lets filtering_r.pl take a variant directly (-v) instead of a pre-made
# annotated VCF. Coordinates are normalized offline; HGVS is recoded to genomic
# coordinates via the Ensembl REST API (only the variant string is transmitted,
# never patient data). The resulting *.germline.vep.vcf.gz is then analyzed in
# lookup mode. Needs curl + jq for the HGVS path (coords need neither).
#############################################################################

# Minimal sites-only VCF header with GRCh38 'chr' contigs (so bcftools norm in
# vep_annotate.sh has them defined).
sub lookup_vcf_header {
    return <<'EOF';
##fileformat=VCFv4.2
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
}

# Resolve one variant string to (chr, pos, REF, ALT). Dashed/colon/space coords
# (chr2-166199981-A-G | 2:166199981:A:G | 'chr2 166199981 A G') are parsed
# offline; anything else is treated as HGVS and recoded via Ensembl REST.
sub resolve_variant {
    my ($v, $work) = @_;
    (my $norm = $v) =~ s/[:\-]/ /g;
    my @f = split ' ', $norm;
    # chr-pos-ref-alt (4 fields) — VCF-style, POS = first REF base.
    if (@f == 4 && $f[1] =~ /^\d+$/
        && $f[2] =~ /^[ACGTNacgtn]+$/ && $f[3] =~ /^[ACGTNacgtn*]+$/) {
        my $chrom = ($f[0] =~ /^chr/i) ? $f[0] : "chr$f[0]";
        printf "[lookup] coordinates: %s:%s %s>%s\n", $chrom, $f[1], uc $f[2], uc $f[3];
        return ($chrom, $f[1], uc $f[2], uc $f[3]);
    }
    # chr-start-end-ref-alt (5 fields) — .candidatos column order; POS = start,
    # END ignored (it's start+len(REF)-1). vep_annotate.sh left-aligns with
    # `bcftools norm -f`, so a VCF-style anchored REF/ALT is what's expected here.
    if (@f == 5 && $f[1] =~ /^\d+$/ && $f[2] =~ /^\d+$/
        && $f[3] =~ /^[ACGTNacgtn]+$/ && $f[4] =~ /^[ACGTNacgtn*]+$/) {
        my $chrom = ($f[0] =~ /^chr/i) ? $f[0] : "chr$f[0]";
        printf "[lookup] coordinates: %s:%s-%s %s>%s (POS=%s, END dropped)\n",
               $chrom, $f[1], $f[2], uc $f[3], uc $f[4], $f[1];
        return ($chrom, $f[1], uc $f[3], uc $f[4]);
    }
    return resolve_hgvs($v, $work);
}

# Recode an HGVS string to GRCh38 genomic coords via Ensembl variant_recoder.
# Needs curl + jq. Dies with an actionable message on failure.
sub resolve_hgvs {
    my ($v, $work) = @_;
    print "[lookup] HGVS input: $v\n";
    print "[lookup] resolving via Ensembl REST ($ENSEMBL_REST) — only the variant notation is sent, no patient data.\n";
    (my $enc = $v) =~ s/([^A-Za-z0-9._~-])/sprintf("%%%02X", ord($1))/ge;
    my $url  = "$ENSEMBL_REST/variant_recoder/human/$enc?content-type=application/json&vcf_string=1";
    my $json = "$work/rest.json";
    system("timeout 30 curl -sS '$url' > '$json' 2>'$work/rest.err'");
    my $resp = "";
    if (open my $jh, "<", $json) { local $/; $resp = <$jh>; close $jh; }
    if (!defined $resp || $resp eq "") {
        warn "[lookup] no response from Ensembl REST (offline, or endpoint unreachable).\n";
        warn "         On an air-gapped host, pass genomic coordinates (chr-pos-ref-alt) instead.\n";
        die  "[lookup] HGVS resolution failed for: $v\n";
    }
    my $err = qx(jq -r 'if type=="object" and has("error") then .error else empty end' '$json' 2>/dev/null);
    chomp $err;
    if ($err ne "") {
        warn "[lookup] Ensembl could not resolve this HGVS: $err\n";
        warn "         Use a transcript-qualified ENST\x{2026} HGVS (this cache is Ensembl, not RefSeq), or pass coordinates.\n";
        die  "[lookup] HGVS resolution failed for: $v\n";
    }
    my $vcfstr = qx(jq -r '.[] | to_entries[] | .value.vcf_string[]?' '$json' 2>/dev/null | grep -E '^(MT|X|Y|[0-9]+)-' | head -1);
    chomp $vcfstr;
    die "[lookup] could not extract genomic coordinates from the HGVS resolution for: $v\n" if $vcfstr eq "";
    my ($rchr,$rpos,$rref,$ralt) = split /-/, $vcfstr;
    $rchr = "M" if $rchr eq "MT";
    $rchr = "chr$rchr";
    printf "[lookup] resolved to: %s:%s %s>%s\n", $rchr, $rpos, $rref, $ralt;
    return ($rchr, $rpos, $rref, $ralt);
}

# Build a sites-only VCF for all requested variants and annotate it through
# vep_annotate.sh. Returns (annotated_vcf, scratch_dir). The annotated VCF is
# named Lookup.<tag>.germline.vep.vcf.gz so lookup mode derives a clean output
# name (Lookup.<tag>.<panel>.candidatos).
sub build_and_annotate_lookup {
    my ($vars) = @_;
    my @reqs = grep { defined && $_ ne "" } @$vars;
    die "no variants given (use -v/--variant)\n" unless @reqs;

    my $work = qx(mktemp -d); chomp $work;
    die "could not create scratch dir (mktemp)\n" unless $work ne "" && -d $work;
    my $raw = "$work/variant.vcf";
    open my $rfh, ">", $raw or die "$raw: $!";
    print $rfh lookup_vcf_header();
    my @ids;
    for my $v (@reqs) {
        my ($chr,$pos,$ref,$alt) = resolve_variant($v, $work);
        printf $rfh "%s\t%s\t.\t%s\t%s\t.\t.\t.\n", $chr, $pos, $ref, $alt;
        push @ids, "$chr-$pos-$ref-$alt";
    }
    close $rfh;

    my $tag = (@ids == 1) ? $ids[0] : $ids[0] . "+" . (scalar(@ids) - 1);
    $tag =~ s/[^A-Za-z0-9._-]/_/g;
    my $ann = "Lookup.$tag.germline.vep.vcf.gz";

    printf "[lookup] %d variant(s) -> sites-only VCF; annotating via vep_annotate.sh (this can take a few minutes) ...\n", scalar @ids;
    my $rc = system("bash '$CF_REPO/vep_annotate.sh' '$raw' '$ann' > '$work/annotate.log' 2>&1");
    if ($rc != 0) {
        system("tail -25 '$work/annotate.log' 1>&2");
        die "[lookup] annotation failed (see vep_annotate.sh output above)\n";
    }

    # Splice scoring: a single-variant consult is meant to report EVERYTHING, so
    # (unless --no-splice) run Pangolin on the variant(s) inline and drop the
    # scores where the main loop expects them (Lookup.<tag>.<panel>.pangolin.tsv).
    # Degrades gracefully — a missing env/reference or a Pangolin failure just
    # leaves pangolin_score blank; it never aborts the consult.
    unless ($NO_SPLICE) {
        my $tsv = "Lookup.$tag.$PANEL_TAG.pangolin.tsv";
        run_pangolin_lookup($ann, $work, $tsv);
    }
    return ($ann, $work);
}

# Score the annotated lookup variant(s) with Pangolin and write parse_pangolin's
# per-variant max|delta| map to $tsv (keyed chr-pos-ref-alt, matching the main
# loop's $my_id — built from the NORMALIZED annotated VCF, not the raw input).
# Returns 1 on success; on any problem warns and returns 0 (pangolin_score stays
# blank). Env overrides mirror run_filtering.sh: $CONDA_BASE / $PANGOLIN_ENV /
# $PANGOLIN_FASTA / $PANGOLIN_DB.
sub run_pangolin_lookup {
    my ($ann, $work, $tsv) = @_;
    my $conda = $ENV{CONDA_BASE}     // "$ENV{HOME}/miniconda3";
    my $env   = $ENV{PANGOLIN_ENV}   // "pangolin";
    my $fa    = $ENV{PANGOLIN_FASTA} // "$ENV{HOME}/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa";
    my $db    = $ENV{PANGOLIN_DB}    // "$ENV{HOME}/vep_refs/pangolin/gencode.v38.annotation.db";
    unless (-e "$conda/etc/profile.d/conda.sh" && -e $fa && -e $db) {
        warn "[lookup] Pangolin skipped — conda env or reference not found (pangolin_score stays blank; use --no-splice to silence).\n";
        return 0;
    }

    my $csv = "$work/pangolin_input.csv";
    open my $cf, ">", $csv or do { warn "[lookup] Pangolin: cannot write $csv ($!)\n"; return 0; };
    print $cf "CHROM,POS,REF,ALT\n";
    my $n = 0;
    my $fh = open_vcf($ann);
    while (my $line = <$fh>) {
        next if $line =~ /^#/;
        my @c = split /\t/, $line;
        next if $c[4] =~ /,/;                       # pre-split; skip any residual multiallelic
        print $cf "$c[0],$c[1],$c[3],$c[4]\n"; $n++;
    }
    close $fh; close $cf;
    return 0 unless $n;

    print "[lookup] scoring splicing with Pangolin ($n variant(s), env '$env') ...\n";
    my $out = "$work/pangolin_out";                  # pangolin writes <out>.csv
    my $sh  = "source '$conda/etc/profile.d/conda.sh' && conda activate '$env' && "
            . "pangolin '$csv' '$fa' '$db' '$out' -c CHROM,POS,REF,ALT";
    if (system("bash", "-c", "$sh > '$work/pangolin.log' 2>&1") != 0 || !-e "$out.csv") {
        warn "[lookup] Pangolin failed (pangolin_score stays blank); tail of log:\n";
        system("tail -12 '$work/pangolin.log' 1>&2");
        return 0;
    }
    if (system("perl '$CF_REPO/parse_pangolin.pl' '$out.csv' > '$tsv' 2>'$work/parse.log'") != 0) {
        warn "[lookup] parse_pangolin.pl failed (pangolin_score stays blank)\n";
        unlink $tsv;
        return 0;
    }
    return 1;
}

#############################################################################
# Discover trios / duos from filenames
#
# Role-suffix naming convention (single source of truth — see sample_role):
#     <FAMILY>-P = proband, <FAMILY>-M = mother, <FAMILY>-F = father.
# Each role shares the FAMILY prefix, so <FAMILY>-P/-M/-F form one trio. A
# sample whose name does not end in -P/-M/-F is ignored by auto-discovery (it
# can still be analyzed via --proband). See run_naming_selftest for examples.
#############################################################################

# Parse a sample base-name into (role, family). Role is 'P' (proband),
# 'M' (mother), 'F' (father), or '' for a non-conforming name.
sub sample_role {
    my ($s) = @_;
    return ($2, $1) if defined $s && $s =~ /^(.+)-([PMF])$/;
    return ("", "");
}

# Group sample base-names into families by the -P/-M/-F convention. Returns an
# arrayref of {proband, mother, father} records (mother/father undef if absent),
# one per family that has a proband, ordered by family name.
sub discover_families {
    my %present = map { $_ => 1 } @_;
    my %fam;
    for my $s (keys %present) {
        my ($role, $f) = sample_role($s);
        $fam{$f}{$role} = $s if $role;
    }
    my @recs;
    for my $f (sort keys %fam) {
        next unless $fam{$f}{P};
        push @recs, { proband => $fam{$f}{P}, mother => $fam{$f}{M}, father => $fam{$f}{F} };
    }
    return \@recs;
}

# Self-test of the naming logic above — no VCFs or reference files required.
sub run_naming_selftest {
    my @ok;
    my $is = sub {
        my ($got, $want, $msg) = @_;
        my ($g, $w) = (defined $got ? $got : "<undef>", defined $want ? $want : "<undef>");
        my $pass = ($g eq $w);
        printf "  [%s] %-26s got=%s want=%s\n", $pass ? "PASS" : "FAIL", $msg, $g, $w;
        push @ok, $pass;
    };

    print "naming self-test\n";
    my @sr;
    @sr = sample_role("EPID107-P"); $is->("$sr[0]/$sr[1]", "P/EPID107", "sample_role proband");
    @sr = sample_role("EPID107-M"); $is->("$sr[0]/$sr[1]", "M/EPID107", "sample_role mother");
    @sr = sample_role("EPID107-F"); $is->("$sr[0]/$sr[1]", "F/EPID107", "sample_role father");
    @sr = sample_role("EPID107");   $is->("$sr[0]/$sr[1]", "/",         "non-conforming (no role)");

    my $recs = discover_families(qw(
        EPID107-P EPID107-M EPID107-F   EPIC280-P EPIC280-M   STRAY junk-X
    ));
    my %by = map { $_->{proband} => $_ } @$recs;
    $is->(scalar @$recs,            2,            "two probands discovered");
    $is->($by{"EPID107-P"}{mother}, "EPID107-M",  "trio mother");
    $is->($by{"EPID107-P"}{father}, "EPID107-F",  "trio father");
    $is->($by{"EPIC280-P"}{mother}, "EPIC280-M",  "duo mother");
    $is->($by{"EPIC280-P"}{father}, undef,        "duo has no father");

    my $fail = grep { !$_ } @ok;
    print $fail ? "naming self-test: $fail FAILED\n" : "naming self-test: all ".scalar(@ok)." passed\n";
    exit($fail ? 1 : 0);
}

# Self-test of the cohort recurrent-artifact filter [#11] — builds a synthetic
# 12-sample cohort of tiny sites+GT VCFs (no VEP/reference files needed), then
# checks the carrier tally and the drop decision, including the founder-safety
# guard (recurrent-but-in-gnomAD is KEPT) and the small-cohort guard (N<min = off).
sub run_cohort_selftest {
    print "cohort recurrent-artifact self-test\n";
    my @ok;
    my $is = sub {
        my ($got,$want,$msg) = @_;
        my $pass = ($got eq $want);
        printf "  [%s] %-38s got=%s want=%s\n", $pass ? "PASS" : "FAIL", $msg, $got, $want;
        push @ok, $pass;
    };

    my $dir = qx(mktemp -d); chomp $dir;
    die "cohort self-test: mktemp failed\n" unless $dir ne "" && -d $dir;
    my $hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
    # 12-sample cohort:
    #   ART_HOM  chr1-100-A-G : hom  in samples 1..6  (6/12 = 50%)  -> recurrent artifact
    #   ART_HET  chr2-200-C-T : het  in samples 1..12 (12/12 = 100%) -> recurrent artifact
    #   REAL     chr3-300-G-A : het  in samples 1..2  (2/12 = 17%)  -> not recurrent
    for my $s (1..12) {
        open my $vh, ">", "$dir/S$s.vcf" or die "cohort self-test: $!";
        print $vh $hdr;
        print $vh "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t1/1\n" if $s <= 6;
        print $vh "chr2\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1\n";
        print $vh "chr3\t300\t.\tG\tA\t.\t.\t.\tGT\t0/1\n" if $s <= 2;
        close $vh;
    }
    my @files = sort glob("$dir/*.vcf");
    my ($ac,$hom,$het,$n) = build_cohort_tally(\@files);
    $is->($n,                            12, "cohort size N");
    $is->($ac->{"chr1-100-A-G"}  // 0,    6, "ART_HOM carriers");
    $is->($hom->{"chr1-100-A-G"} // 0,    6, "ART_HOM all homozygous");
    $is->($het->{"chr1-100-A-G"} // 0,    0, "ART_HOM no hets");
    $is->($ac->{"chr2-200-C-T"}  // 0,   12, "ART_HET carriers");
    $is->($ac->{"chr3-300-G-A"}  // 0,    2, "REAL carriers");
    # Drop decision (3rd arg is the gnomAD joint AC; absence means AC == 0).
    $is->(cohort_artifact_call($ac->{"chr1-100-A-G"}, $n, 0),    1, "ART_HOM dropped (recurrent+absent)");
    $is->(cohort_artifact_call($ac->{"chr2-200-C-T"}, $n, 0),    1, "ART_HET dropped (recurrent+absent)");
    $is->(cohort_artifact_call($ac->{"chr3-300-G-A"}, $n, 0),    0, "REAL kept (not recurrent)");
    $is->(cohort_artifact_call($ac->{"chr1-100-A-G"}, $n, 1),    0, "recurrent but AC=1 in gnomAD -> kept (founder-safe)");
    $is->(cohort_artifact_call($ac->{"chr1-100-A-G"}, $n, 800),  0, "recurrent but common in gnomAD -> kept");
    $is->(cohort_artifact_call(3, 4, 0),                         0, "small cohort (N<min) -> filter off");
    # The absolute carrier floor is what protects small batches: at N=8 the 25% fraction
    # alone would drop a 2-carrier variant, which is an ordinary shared/founder allele.
    $is->(cohort_artifact_call(2, 8, 0),                         0, "2 of 8 carriers -> below floor, kept");
    $is->(cohort_artifact_call(3, 8, 0),                         1, "3 of 8 carriers -> at floor, dropped");
    $is->(cohort_artifact_call(8, 8, 0),                         1, "8 of 8 carriers -> dropped");
    # ...and the fraction is what protects large cohorts, where 3 carriers is noise.
    $is->(cohort_artifact_call(3, 56, 0),                        0, "3 of 56 carriers -> below fraction, kept");
    $is->(cohort_artifact_call(14, 56, 0),                       1, "14 of 56 carriers -> dropped");
    system("rm","-rf",$dir);

    my $fail = grep { !$_ } @ok;
    print $fail ? "cohort self-test: $fail FAILED\n" : "cohort self-test: all ".scalar(@ok)." passed\n";
    exit($fail ? 1 : 0);
}

my @files = glob $INPUT_GLOB;
my %file_for;
for my $f (@files) { (my $s = $f) =~ s/\.germline\.vep\.vcf\.gz$//; $file_for{$s} = $f; }

# Proband list + parent map: forced override (--proband) or auto-discovery.
my @probands;
my %parents_of;        # proband -> {m=>base|undef, f=>base|undef}
if (@force_probands) {
    for my $p (@force_probands) {
        die "forced proband '$p' has no '$p.germline.vep.vcf.gz' in this directory\n"
            unless exists $file_for{$p};
        (my $fam = $p) =~ s/-P$//;
        $parents_of{$p} = { m => (exists $file_for{"$fam-M"} ? "$fam-M" : undef),
                            f => (exists $file_for{"$fam-F"} ? "$fam-F" : undef) };
    }
    @probands = @force_probands;
} else {
    for my $r (@{ discover_families(keys %file_for) }) {
        push @probands, $r->{proband};
        $parents_of{$r->{proband}} = { m => $r->{mother}, f => $r->{father} };
    }
    # Singleton fallback: auto-discovery keys on the -P/-M/-F role suffix. If NO family
    # structure is found but samples exist, treat every sample as a standalone proband
    # (parents still auto-linked if -M/-F siblings happen to exist). Prevents a SILENT
    # 0-candidatos run on cohorts named without role suffixes (e.g. EPIGEN01..20). Only
    # fires when discovery found nothing, so trio/duo runs are unaffected.
    if (!@probands && keys %file_for && !$LOOKUP && !@VARIANTS) {
        warn "NOTE: no -P/-M/-F family structure found; analyzing all "
           . (scalar keys %file_for) . " sample(s) as singleton probands "
           . "(inheritance=NA; use -P/-M/-F names or --proband for trio/duo analysis).\n";
        for my $p (sort keys %file_for) {
            push @probands, $p;
            (my $fam = $p) =~ s/-P$//;
            $parents_of{$p} = { m => (exists $file_for{"$fam-M"} ? "$fam-M" : undef),
                                f => (exists $file_for{"$fam-F"} ? "$fam-F" : undef) };
        }
    }
}

# ── Variant entry point (-v): build + annotate, then enter lookup mode ──
# Resolve the requested variant(s) to a sites-only chr-prefixed VCF, annotate it
# through vep_annotate.sh, and feed the result into lookup mode below. The
# annotated VCF + scratch dir are cleaned up after the run (unless --keep-vcf).
my ($LOOKUP_ANN, $LOOKUP_WORK);
if (@VARIANTS) {
    die "give EITHER variant input (-v) OR --lookup, not both\n" if $LOOKUP;
    ($LOOKUP_ANN, $LOOKUP_WORK) = build_and_annotate_lookup(\@VARIANTS);
    $LOOKUP      = 1;
    $LOOKUP_FILE = $LOOKUP_ANN;
}

# ── Single-variant lookup mode ──
# One explicit annotated VCF (sites-only or single-sample), analyzed as a
# singleton with the panel / rarity / consequence / rescue gates and the
# recessive-carrier drop all bypassed: report EVERY qualifying annotation of
# the variant in the transposed readable format. MANE-only unless
# --all-transcripts. Overrides any filename-based discovery above.
if ($LOOKUP) {
    die "--lookup needs an annotated VCF (<name>.germline.vep.vcf.gz)\n"
        unless defined $LOOKUP_FILE && $LOOKUP_FILE ne "";
    die "--lookup file not found: $LOOKUP_FILE\n" unless -e $LOOKUP_FILE;
    (my $name = $LOOKUP_FILE) =~ s{.*/}{};
    $name =~ s/\.germline\.vep\.vcf\.gz$//;
    $name =~ s/\.vcf(\.gz)?$//;
    %file_for   = ($name => $LOOKUP_FILE);
    @probands   = ($name);
    %parents_of = ($name => { m => undef, f => undef });
}

print "muestras encontradas: ", join(", ", sort keys %file_for), "\n";
print "probandos: ", join(", ", @probands), (@force_probands ? " (forced override)" : ""), "\n";

# Output column order (single definition, reused for header + rows).
my @COLS = qw(
    chr start end ref alt gene strand consequence hgvs
    revel eve_class eve_score cadd am_class am_score pangolin_score
    clinvar_sig clinvar_stars clinvar_disease clinvar_aa
    loftee
    gnomAD_ac gnomAD_an gnomAD_af gnomAD_nhomalt gnomAD_filter
    zygosity GT DP GQ AB GT_SOURCE NCALLERS CONF
    inheritance kept_by acmg_class acmg_criteria flags
    Association MOI GDV
);

# Cohort recurrent-artifact state [#11] — built lazily on the first FINAL-pass
# proband (see below), and only for a real cohort auto-analyzed together. A single
# variant (-v/--lookup), a forced/single proband (--proband), or a run of fewer than
# $COHORT_MIN probands has no cohort to compare against, so the filter stays OFF.
my ($cohort_ac,$cohort_hom,$cohort_het,$cohort_n) = ({},{},{},0);
my $cohort_on    = 0;
my $cohort_built = 0;

# Every final-pass candidate row from every proband, each prefixed with its sample id.
# Written once, after the loop, as batch.<panel>.candidatos — one file to open when the
# question is "what did this batch turn up" rather than "what did this patient turn up".
my @batch_rows;

#############################################################################
# Process each proband
#############################################################################

foreach my $proband (@probands) {
    my $pfile = $file_for{$proband};
    my $mname = $parents_of{$proband}{m};
    my $fname = $parents_of{$proband}{f};
    my $mfile = defined $mname ? $file_for{$mname} : undef;
    my $ffile = defined $fname ? $file_for{$fname} : undef;
    my $have_m = defined $mfile;
    my $have_f = defined $ffile;

    my $tsv   = "$proband.$PANEL_TAG.pangolin.tsv";
    my $final = $LOOKUP ? 1 : -e $tsv;          # lookup reports immediately (no 2-pass)
    my $pscore = (-e $tsv) ? load_scores($tsv) : {};

    # ── [#11] Build the cohort recurrent-artifact tally once, on the first FINAL-pass
    #    proband. Gated to a real cohort auto-analyzed together: OFF for lookup, for a
    #    forced/single proband (--proband), and for runs of < $COHORT_MIN probands. ──
    if ($final && !$LOOKUP && !$cohort_built) {
        $cohort_built = 1;
        # Eligibility is counted in SAMPLES, the same unit as the denominator, so the
        # acting threshold matches the one documented. (Counting probands for
        # eligibility while dividing by all VCFs silently diluted the fraction ~3x on
        # trio cohorts: an artifact in 4 of 10 probands scored 4/30 = 13%, under the
        # 25% bar.) Parents belong in the tally — a mapping artifact is a property of
        # the assay, not of affected status, so every extra sample sharpens it.
        my $eligible = (!@force_probands && @probands >= 2 && @files >= $COHORT_MIN);
        if ($eligible) {
            ($cohort_ac,$cohort_hom,$cohort_het,$cohort_n) = build_cohort_tally(\@files);
            $cohort_on = ($cohort_n >= $COHORT_MIN) ? 1 : 0;
        }
        if ($cohort_on) {
            printf "cohort artifact filter: ON — N=%d samples; drop if carried by >=%d samples AND >=%d%% of cohort AND absent from gnomAD (AC=0)%s\n",
                   $cohort_n, $COHORT_MIN_CARRIERS, int($COHORT_MAX_FRAC*100 + 0.5),
                   $KEEP_COHORT_ARTIFACTS ? " (keep+tag mode)" : "";
        } else {
            printf "cohort artifact filter: OFF — %s\n",
                   (@force_probands ? "forced/single proband (--proband)"
                    : @probands < 2 ? "single proband"
                    : "only ".scalar(@files)." sample(s); need >=$COHORT_MIN");
        }
    }

    my $kind = ($have_m && $have_f) ? "trio" : ($have_m||$have_f) ? "duo" : "singleton";
    my $final_desc = $final ? ((-e $tsv) ? "FINAL: splice scores from $tsv" : "FINAL: no splice scores")
                            : "EMIT: no scores yet";
    print "\ntrabajando con $proband ($kind, $final_desc)";
    print " | madre: $mfile" if $have_m;
    print " | padre: $ffile" if $have_f;
    print "\n";

    # Resolve CSQ indices for this proband's layout.
    my $col = csq_columns($pfile);
    my %i = (
        gene          => resolve($col,'SYMBOL'),
        strand        => resolve($col,'STRAND'),
        transcript    => resolve($col,'Feature'),
        mane_select   => resolve($col,'MANE_SELECT'),
        consequence   => resolve($col,'Consequence'),
        hgvsc         => resolve($col,'HGVSc'),
        hgvsp         => resolve($col,'HGVSp'),
        tpos          => resolve($col,'cDNA_position'),
        aa            => resolve($col,'Amino_acids'),
        ppos          => resolve($col,'Protein_position'),
        revel         => resolve($col,'REVEL'),
        eve_class     => resolve($col,'EVE_CLASS'),
        eve_score     => resolve($col,'EVE_SCORE'),
        cadd          => resolve($col,'CADD_PHRED'),
        am_class      => resolve($col,'am_class'),
        am_score      => resolve($col,'am_pathogenicity'),
        loftee        => resolve($col,'LoF'),
        loftee_filter => resolve($col,'LoF_filter'),
        loftee_flags  => resolve($col,'LoF_flags'),
        g_ac          => resolve($col,'gnomADmin_AC_joint','gnomad.*AC'),
        g_an          => resolve($col,'gnomADmin_AN_joint','gnomad.*AN'),
        g_nhom        => resolve($col,'gnomADmin_nhomalt_joint','gnomad.*nhomalt'),
        g_filter      => resolve($col,'gnomADmin_FILTER','gnomad.*FILTER'),
        # Fresh ClinVar (--custom) preferred; fall back to cache CLIN_SIG.
        clnsig        => resolve($col,'ClinVar_CLNSIG','CLIN_SIG'),
        clnstars      => resolve($col,'ClinVar_CLNREVSTAT'),
        clndn         => resolve($col,'ClinVar_CLNDN'),
        # PERv1 pathogenic-variant-enriched region overlap, for ACMG PM1. Optional:
        # absent -> PM1 never fires (no warning; not every run annotates it).
        per           => resolve($col,'PER'),
    );

    # [#3] Assert critical fields resolved — fail loudly, never silently pass-all.
    my @critical = qw(gene transcript consequence cadd g_ac g_an);
    my @missing = grep { !defined $i{$_} } @critical;
    die "FATAL: CSQ fields not found in $pfile: @missing\n".
        "  The annotation layout changed; check the CSQ header.\n" if @missing;
    warn "WARN: ClinVar field not found in $pfile (ClinVar rescue/columns disabled)\n"
        unless defined $i{clnsig};

    my ($mama,$papa) = ({},{});
    if ($final) { $mama = load_parent($mfile); $papa = load_parent($ffile); }

    # Run statistics [#9].
    my %stat = (lines=>0, multiallelic=>0, structural=>0, cohort_dropped=>0, probes=>0);

    my %emit;            # EMIT pass: unique candidate variants
    my %emit_probe;      # vid -> 1 for splice-discovery probes (counted once) [#5]
    my @rows;            # FINAL pass: buffered rows (for per-gene recessive logic)
    my %cohort_seen;     # [#11] vid -> 1 once dropped/logged (avoid double-counting)

    my $pfh = open_vcf($pfile);
    while (my $v = <$pfh>) {
        chomp $v;
        next if $v =~ /^#/;
        $stat{lines}++;

        my @c = split /\t/, $v;
        my ($chr,$start,$ref,$alt,$info,$fmt,$smp) = @c[0,1,3,4,7,8,9];

        if ($alt =~ /,/) { $stat{multiallelic}++; next; }   # [#2] should be pre-split

        my ($csq) = $info =~ /(?:^|;)CSQ=([^;]*)/;
        # PERv1 overlap may also sit in a plain INFO tag when the BED was applied
        # positionally (outside VEP). Read it here so the per-transcript loop below
        # can fall back to it; the region still has to name this record's gene.
        my ($per_info) = $info =~ /(?:^|;)PER=([^;]*)/;
        next unless defined $csq;

        my $my_id = "$chr-$start-$ref-$alt";

        # Consensus provenance (INFO-level; present only for consensus.sh output —
        # DRAGEN VCFs lack these, leaving the columns empty). GT_SOURCE flags
        # genotypes borrowed from a non-DeepVariant caller (which also lack VAF).
        my ($gtsrc)    = $info =~ /(?:^|;)GT_SOURCE=([^;]*)/; $gtsrc    = defined $gtsrc    ? $gtsrc    : "";
        my ($ncallers) = $info =~ /(?:^|;)NCALLERS=([^;]*)/;  $ncallers = defined $ncallers ? $ncallers : "";
        my ($conf)     = $info =~ /(?:^|;)CONF=([^;]*)/;      $conf     = defined $conf     ? $conf     : "";

        # Proband genotype (same for all transcripts of this variant) [#5].
        my ($gt,$dp,$gq,$adr,$ada) = parse_call($fmt,$smp);
        my $zyg = zygosity($gt);
        my $ab  = ($ada ne "" && ($adr+$ada) > 0) ? sprintf("%.2f", $ada/($adr+$ada)) : "";

        foreach my $fila (split /,/, $csq) {
            my @r = split /\|/, $fila, -1;
            my $gene        = field(\@r,$i{gene});
            my $transcript  = field(\@r,$i{transcript});
            my $consequence = field(\@r,$i{consequence});

            next if !exists($mane{$transcript}) && !($LOOKUP && $ALL_TX);  # MANE-only unless --all-transcripts (lookup)
            my $in_panel = exists $epigenes{$gene};
            my $in_acmg  = exists $acmg{$gene};
            next if !($in_panel || $in_acmg) && !$LOOKUP;   # lookup reports off-panel genes too

            my $g_ac = field(\@r,$i{g_ac}); $g_ac = ($g_ac eq "") ? 0 : $g_ac;
            my $g_an = field(\@r,$i{g_an}); $g_an = ($g_an eq "") ? 0 : $g_an;
            my $freq = ($g_an > 0) ? ($g_ac/$g_an)*100 : 0;

            my $clnsig    = field(\@r,$i{clnsig});
            my $clnstars  = field(\@r,$i{clnstars});
            my $clndn     = field(\@r,$i{clndn});

            # An established ClinVar classification (P/LP with >=1 review star) is
            # itself evidence, so it admits a panel-gene variant regardless of the MOI
            # frequency ceiling and the consequence whitelist. Both gates otherwise
            # fire BEFORE any evidence arm is consulted, silently discarding two
            # classes of true positive: founder alleles that sit above the ceiling in
            # a bottlenecked population, and pathogenic non-coding / synonymous
            # variants whose consequence term is not whitelisted. The ACMG-SF path
            # exempts the same tier for the same reason.
            my $clinvar_established =
                (clinvar_pathogenic($clnsig) && clinvar_stars($clnstars) >= 1) ? 1 : 0;

            # HGVSc is needed BEFORE the gate (the splice probe reads its intron offset),
            # so it is extracted here rather than in the final-pass block below.
            my $hgvsc = field(\@r,$i{hgvsc});

            # Candidate structural gate: panel gene + whitelisted consequence + rare
            # (MOI-aware: recessive genes tolerate higher carrier freq). [#1,#6]
            my $cand_structural = 0;
            my ($p_assoc,$p_moi,$p_gdv) = ("","","");
            if ($in_panel) { ($p_assoc,$p_moi,$p_gdv) = split /\t/, $epigenes{$gene}; $p_moi //= ""; }
            if ($in_panel && (consequence_ok($consequence) || $clinvar_established)) {
                my $recessive = moi_recessive($p_moi);
                $cand_structural = ($clinvar_established
                                 || $freq <= ($recessive ? $FREQ_AR : $FREQ_AD)) ? 1 : 0;
            }
            if ($LOOKUP) { $cand_structural = 1; }   # report-everything: force structural pass

            # [#6] ACMG-SF genes must reach Pangolin too. The SF tier itself doesn't use a
            # splice score, but without this every GDV=Incidental row carried a blank
            # pangolin_score and could never earn BP7 or the splice rescue — a secondary
            # finding was structurally denied the evidence a primary candidate gets.
            my $sf_structural = 0;
            if ($in_acmg && !$cand_structural
                && (consequence_ok($consequence) || $clinvar_established)) {
                $sf_structural = ($clinvar_established || $freq <= $SF_FREQ_MAX) ? 1 : 0;
            }

            # [#5] Splice DISCOVERY probe — an intronic/synonymous variant that no gate
            # would ever admit, scored so Pangolin can find a splice disruption in it.
            # Bounded by intron distance and a strict rarity ceiling (see $SPLICE_PROBE).
            # $g_an > 0 means the position actually EXISTS in the gnomAD resource. This is
            # load-bearing, not a nicety. The custom VCF is gnomAD.joint.v4.1.**mane**, i.e.
            # MANE transcripts plus flanks — deep intronic sequence is largely absent from
            # it. A variant there gets AC=""/AN="" which this code coerces to 0, so $freq
            # computes as 0 and passes ANY rarity ceiling, and PM2 ("absent or singleton in
            # gnomAD") fires on what is really an annotation gap. Probing uncovered
            # territory would therefore rescue variants with no working frequency filter and
            # a manufactured pathogenic criterion. Requiring coverage bounds the probe set
            # to where the resource can actually answer the question.
            my $probe = 0;
            if ($SPLICE_PROBE && !$LOOKUP && !$cand_structural && !$sf_structural
                && ($in_panel || $in_acmg)
                && (!$PROBE_REQUIRE_GNOMAD || $g_an > 0)
                && $freq <= $PROBE_FREQ_MAX && !consequence_ok($consequence)) {
                if ($consequence =~ /(?:^|&)synonymous_variant(?:&|$)/) {
                    $probe = 1;                       # exonic splice-altering synonymous
                } elsif ($consequence =~ /(?:^|&)intron_variant(?:&|$)/) {
                    my $off = intron_offset($hgvsc);
                    $probe = 1 if defined $off && $off <= $INTRON_MAX_DIST;
                }
            }

            # Collect variants for Pangolin. Done in BOTH passes so the input CSV is
            # (re)written with the full set even in the final pass (otherwise it would be
            # emptied and the score cache wiped).
            $emit{$my_id} = "$chr,$start,$ref,$alt"
                if ($cand_structural || $sf_structural || $probe) && !exists $emit{$my_id};
            $stat{probes}++ if $probe && !$emit_probe{$my_id}++;
            next unless $final;
            $stat{structural}++ if $cand_structural;

            # ── FINAL pass: extract scoring fields (shared by both paths) ──
            my $revel     = field(\@r,$i{revel});
            my $eve_class = field(\@r,$i{eve_class});
            my $eve_score = field(\@r,$i{eve_score});
            my $cadd      = field(\@r,$i{cadd});
            my $am_class  = field(\@r,$i{am_class});
            my $am_score  = field(\@r,$i{am_score});
            my $loftee    = field(\@r,$i{loftee});
            my $lof_filt  = field(\@r,$i{loftee_filter});
            my $lof_flag  = field(\@r,$i{loftee_flags});
            my $g_nhom    = field(\@r,$i{g_nhom});
            my $g_filter  = field(\@r,$i{g_filter});
            my $strand    = field(\@r,$i{strand});
            my $hgvsp     = field(\@r,$i{hgvsp});
            my $tpos      = field(\@r,$i{tpos});
            my $clnstar_n = clinvar_stars($clnstars);          # exact-variant review stars
            my $pangolin  = exists $pscore->{$my_id} ? $pscore->{$my_id} : "";

            # ── PS1 / PM5 from ClinVar amino-acid evidence (>=1 star) ──
            # PS1: SAME missense change is P/LP. PM5: a DIFFERENT change at the same
            # residue is P/LP. EXTENDED: a single-codon in-frame deletion removing
            # residue X is treated as PM5 when any missense at X is P/LP — a
            # different protein change at the same residue is pathogenic, i.e. the
            # residue is intolerant (per-ACMG this is a curatorial extension of PM5
            # beyond missense; flagged in the detail as "(in-frame del)"). Flagged
            # (conflicting) when the matched change is also B/LB. Needs the AA resource.
            my ($aa_crit,$aa_detail,$aa_conflict) = ("","",0);
            my $aa_field = field(\@r,$i{aa});                  # missense "R/H"; 1-codon del "I/-"
            my $ppos_raw = field(\@r,$i{ppos});                # "1922" (single) or "1922-1924" (range)
            my ($aa_ref,$aa_alt) = split m{/}, $aa_field, 2;
            $aa_ref = defined $aa_ref ? uc $aa_ref : "";
            $aa_alt = defined $aa_alt ? uc $aa_alt : "";
            my ($paapos) = ($ppos_raw =~ /(\d+)/);
            $paapos = defined $paapos ? $paapos : "";
            my $single_res  = ($ppos_raw =~ /^\d+$/) ? 1 : 0;                       # one residue, no range
            my $is_missense = ($aa_ref =~ /^[A-Z]$/ && $aa_alt =~ /^[A-Z]$/) ? 1 : 0;
            my $is_codondel = (!$is_missense && $single_res && $aa_ref =~ /^[A-Z]$/
                               && $aa_alt =~ /^-?$/ && $consequence =~ /inframe_deletion/) ? 1 : 0;
            if ($CLINVAR_AA_ON && $gene ne "" && $paapos ne "" && ($is_missense || $is_codondel)) {
                my $refAA = $aa_ref;
                my $altAA = $aa_alt;                       # "" for a deleted residue
                my $k   = "$gene\t$paapos\t$refAA";
                my $plp = $PLP_resid->{$k};
                if ($plp) {
                    # PS1 = a DIFFERENT variant producing the same amino-acid change is
                    # P/LP, so this variant's own ClinVar record is excluded from the
                    # match; self-matching would double-count one submission as PS1+PP5.
                    my $ps1_st = $is_missense ? aa_best_stars($plp->{$altAA}, $my_id) : 0;
                    if ($ps1_st >= 1) {                                                  # PS1
                        $aa_crit   = "PS1";
                        $aa_detail = sprintf("PS1:%s p.%s%s%s (%d*)", $gene,$refAA,$paapos,$altAA,$ps1_st);
                        $aa_conflict = 1 if exists $BLB_resid->{$k}{$altAA};
                    } else {                                                             # PM5
                        my ($balt,$bst) = ("",0);
                        for my $a (keys %$plp) {
                            next if $is_missense && $a eq $altAA;   # exact match is PS1, handled above
                            # A different AA change is necessarily a different variant.
                            my $st = aa_best_stars($plp->{$a}, "");
                            ($balt,$bst) = ($a,$st) if $st >= 1 && $st > $bst;
                        }
                        if ($balt ne "") {
                            $aa_crit = "PM5";
                            $aa_detail = $is_codondel
                              ? sprintf("PM5:%s p.%s%sdel (in-frame del) [P/LP missense at residue: p.%s%s%s (%d*)]",
                                        $gene,$refAA,$paapos, $refAA,$paapos,$balt,$bst)
                              : sprintf("PM5:%s p.%s%s%s [P/LP at residue: p.%s%s%s (%d*)]",
                                        $gene,$refAA,$paapos,$altAA, $refAA,$paapos,$balt,$bst);
                            $aa_conflict = 1 if exists $BLB_resid->{$k}{$balt};
                        }
                    }
                    $aa_detail .= " |conflicting" if $aa_crit && $aa_conflict;
                }
            }
            my $cadd_num  = ($cadd eq "") ? 0 : $cadd;
            my $lof_type  = grep { $LOF_CONS{$_} } split /&/, $consequence;

            # ── Primary candidate inclusion gate (OR) [#4,#8] ──
            my (@kept, $class, $assoc, $moi, $gdv);
            my $sf_ar = 0;
            if ($cand_structural) {
                push @kept, "CADD"     if $cadd_num >= $CADD_MIN;
                push @kept, "AM"       if $am_score ne "" && $am_score >= $AM_MIN;
                push @kept, "EVE"      if $eve_class =~ /athogenic/;
                push @kept, "REVEL"    if $revel ne "" && $revel >= $REVEL_MIN;
                push @kept, "Pangolin" if $pangolin ne "" && $pangolin >= $SPLICE_MIN;
                push @kept, "ClinVar"  if clinvar_pathogenic($clnsig);
                push @kept, "LoF"      if $loftee eq "HC" || ($lof_type && $loftee ne "LC");
                push @kept, $aa_crit   if $aa_crit;   # PS1/PM5 (ClinVar amino-acid evidence)
                # AR_hom: a HOMOZYGOUS rare MANE *protein-altering* variant (missense / inframe /
                # stop_lost / start_lost) in a recessive (AR/XLR) panel gene is rescued even without
                # in-silico/ClinVar support. Rationale (general, not case-specific): a biallelic
                # (homozygous) genotype in a recessive disease gene is itself pathogenicity evidence
                # under recessive inheritance, independent of missense predictors — which are calibrated
                # largely on dominant/heterozygous effects and can miss true recessive alleles. Restricted
                # to coding changes so it does NOT flood on benign homozygous intronic/polypyrimidine
                # variants (those go through the Pangolin splice arm; truncating LoF through the LoF arm);
                # AB>0.75 guards against false-hom calls. Already rare (AR freq gate), MANE, in-panel by
                # this point. BS1/BS2/BA1 flags still annotate benign-leaning ones for the curator.
                # "hem" qualifies on the same rationale: a hemizygous male call in an
                # X-linked recessive gene is the complete genotype, so it is the exact
                # situation the arm exists for. Tagged AR_hem so the curator can see
                # which genotype produced the rescue.
                push @kept, ($zyg eq "hem" ? "AR_hem" : "AR_hom")
                                       if moi_recessive($p_moi) && zyg_biallelic($zyg)
                                          && $consequence =~ /missense|inframe|stop_lost|start_lost|protein_altering/
                                          && $ab ne "" && $ab > 0.75;   # clean hom/hem call (guards against false-hom/artifact; AB~1.0 expected)
                if (@kept) { $class = "primary"; ($assoc,$moi,$gdv) = ($p_assoc,$p_moi,$p_gdv); }
            }

            # ── [#5] Splice-discovery rescue ──
            # A probe is NOT a candidate. It entered the Pangolin set only so the model
            # could look at it, and it earns a row solely on a positive splice score. A
            # probe Pangolin scores low simply vanishes — which is what keeps the probe set
            # from widening the table. Placed before the ACMG-SF block so a splice-active
            # synonymous/intronic variant in a secondary-findings gene is reported as an
            # incidental rather than being tested against the SF coding tiers it cannot meet.
            if (!@kept && $probe && $pangolin ne "" && $pangolin >= $SPLICE_MIN) {
                @kept = ("Pangolin");
                if ($in_panel) {
                    $class = "primary";
                    ($assoc,$moi,$gdv) = ($p_assoc,$p_moi,$p_gdv);
                } else {
                    my ($a_cond,$a_moi) = split /\t/, $acmg{$gene};
                    $class = "incidental";
                    ($assoc,$moi,$gdv) = ($a_cond, $a_moi//"", "Incidental");
                    $sf_ar = (($a_moi//"") =~ /\bAR\b/) ? 1 : 0;
                }
            }

            # ── Incidental (ACMG SF) — stringent; only if not a primary candidate ──
            # ClinVar P/LP (>=1 star) reported regardless of frequency (known founder
            # alleles); novel LoF / >=2-strong-computational tiers require rarity.
            if (!@kept && $in_acmg) {
                my ($a_cond,$a_moi,$a_cat) = split /\t/, $acmg{$gene};
                my $cat_ok =
                    ($a_cat eq "TRUNCATING_ONLY") ? ($lof_type ? 1 : 0)
                  : ($a_cat eq "C282Y_HOM")
                        ? ((($hgvsc =~ /845G>A/ || $hgvsp =~ /Cys282Tyr/i) && $zyg eq "hom") ? 1 : 0)
                  : 1;
                if ($cat_ok) {
                    my @sf;
                    push @sf, "ClinVar_P/LP"
                        if clinvar_pathogenic($clnsig) && clinvar_stars($clnstars) >= 1;
                    if ($freq <= $SF_FREQ_MAX) {
                        push @sf, "LoF" if $loftee eq "HC";
                        my $ncomp = 0;
                        $ncomp++ if $am_score ne "" && $am_score >= $SF_AM;
                        $ncomp++ if $cadd_num >= $SF_CADD;
                        $ncomp++ if $eve_class =~ /athogenic/;
                        $ncomp++ if $revel ne "" && $revel >= $SF_REVEL;
                        push @sf, "Computational" if $ncomp >= 2;
                    }
                    if (@sf) {
                        @kept = @sf; $class = "incidental";
                        ($assoc,$moi,$gdv) = ($a_cond, $a_moi, "Incidental");
                        $sf_ar = ($a_moi =~ /\bAR\b/) ? 1 : 0;
                    }
                }
            }

            # Lookup: report the variant even when no evidence arm fires; supply
            # panel / ACMG-SF / NA context for the Association/MOI/GDV columns.
            if ($LOOKUP && !@kept) {
                $class = "lookup";
                if ($in_panel) {
                    ($assoc,$moi,$gdv) = ($p_assoc,$p_moi,$p_gdv);
                } elsif ($in_acmg) {
                    my ($ac_cond,$ac_moi) = split /\t/, $acmg{$gene};
                    ($assoc,$moi,$gdv) = ($ac_cond, $ac_moi, "Incidental?");
                } else {
                    ($assoc,$moi,$gdv) = ("NA","NA","NA");
                }
            }
            next unless @kept || $LOOKUP;
            my $kept_by = @kept ? join(";", @kept) : "none";

            # Lookup: always surface ACMG-SF context for a secondary-findings gene
            # that isn't in the candidate panel (panel context wins when both).
            if ($LOOKUP && $in_acmg && !$in_panel) {
                my ($ac_cond,$ac_moi) = split /\t/, $acmg{$gene};
                $assoc = $ac_cond     if !defined $assoc || $assoc eq "";
                $moi   = $ac_moi      if !defined $moi   || $moi   eq "";
                $gdv   = "Incidental" if !defined $gdv   || $gdv   eq "";
            }

            # ── Inheritance from parental GT (carrier = non-ref) [#10] ──
            my $in_m = exists $mama->{$my_id};
            my $in_f = exists $papa->{$my_id};
            my $inheritance;
            if ($have_m && $have_f) {
                $inheritance = ($in_m && $in_f) ? "IB" : $in_m ? "IM" : $in_f ? "IF" : "DN";
            } elsif ($have_m) {
                $inheritance = $in_m ? "IM" : "DN/IF";
            } elsif ($have_f) {
                $inheritance = $in_f ? "IF" : "DN/IM";
            } else {
                $inheritance = "NA";
            }

            # ── [#6] Parental call quality for inherited variants ──
            my $inh_lowqual = 0;
            for my $par ($in_m ? $mama->{$my_id} : (), $in_f ? $papa->{$my_id} : ()) {
                my (undef,$pdp,$pgq) = split /:/, $par;
                $inh_lowqual = 1 if ($pdp ne "" && $pdp < $QC_MIN_DP)
                                 || ($pgq ne "" && $pgq < $QC_MIN_GQ);
            }

            # ── [#7] QC / artifact flags ──
            my @qc;
            push @qc, "lowDP"  if $dp ne "" && $dp < $QC_MIN_DP;
            push @qc, "lowGQ"  if $gq ne "" && $gq < $QC_MIN_GQ;
            if ($ab ne "") {
                push @qc, "AB_het" if $zyg eq "het" && ($ab < 0.25 || $ab > 0.75);
                push @qc, "AB_hom" if zyg_biallelic($zyg) && $ab < 0.85;   # hom or hemizygous
            }
            push @qc, "homopolymer"   if homopolymer_context($chr,$start,$ref,$alt);
            push @qc, "GT_rescued"    if $gtsrc ne "" && $gtsrc ne "deepvariant";  # borrowed (non-DV) genotype, no VAF
            push @qc, "inh_lowqual"   if $inh_lowqual;
            push @qc, "DN_unconfirmed" if $inheritance =~ /^DN/;   # no parental ref depth
            my $qc_flag = join(";", @qc);

            # ── [#11] Cohort recurrent-artifact filter (internal panel-of-normals) ──
            # Drop a candidate that is BOTH cohort-recurrent AND gnomAD-absent (a
            # systematic technical artifact; a founder allele would have a gnomAD
            # footprint). --keep-cohort-artifacts keeps it, tagged in `flags` instead.
            if ($cohort_on) {
                my $cac = $cohort_ac->{$my_id} // 0;
                if (cohort_artifact_call($cac, $cohort_n, $g_ac)) {
                    if ($KEEP_COHORT_ARTIFACTS) {
                        push @qc, "cohort_artifact";
                        $qc_flag = join(";", @qc);
                    } else {
                        unless ($cohort_seen{$my_id}++) {
                            $stat{cohort_dropped}++;
                            printf "  cohort_artifact drop: %-28s %-10s carriers=%d/%d (hom=%d het=%d) gnomAD_AC=%s\n",
                                   $my_id, ($gene ne "" ? $gene : "."), $cac, $cohort_n,
                                   ($cohort_hom->{$my_id} // 0), ($cohort_het->{$my_id} // 0),
                                   ($g_ac ne "" ? $g_ac : "0");
                        }
                        next;   # drop this candidate row
                    }
                }
            }

            # ── [#2] Automated ACMG/AMP classification (triage) ──
            # ── PM1: PERv1 region overlap (see acmg_classify) ──
            # VEP emits the BED name field; "/"-delimited, "&"-joined if a variant
            # falls in more than one region. The gene is carried in the name so the
            # region can be checked against the gene of THIS CSQ record: overlapping
            # MANE gene models co-locate in a real exome, and a PER belongs to the
            # gene it was computed on, not to whatever else spans those bases.
            # A positional BED overlap is not transcript-specific, so it may arrive
            # inside CSQ (VEP --custom) or as a plain INFO tag (BED applied outside
            # VEP). MERGE both rather than preferring one: a VCF can carry a partial
            # track in CSQ and the rest in INFO, and preferring CSQ would then hide
            # the INFO regions whenever CSQ happened to be non-empty -- including
            # when its only region belongs to a co-located gene and is rejected by
            # the gene guard below. Duplicates are harmless; ranking dedupes.
            my $per_csq = defined $i{per} ? ($r[$i{per}] // "") : "";
            my $per_raw = join("&", grep { defined && $_ ne "" } ($per_csq, $per_info));
            my ($pm1_strength, $pm1_detail, $pm1_src, $pm1_rank) = ("", "", "", -1);
            if ($per_raw ne "") {
                for my $hit (split /&/, $per_raw) {
                    my @pf = split m{/}, $hit;
                    next unless @pf >= 3;
                    next unless defined $gene && $gene ne "" && $pf[0] eq $gene;
                    my ($st) = $pf[2] =~ /^PM1_(Strong|Moderate)$/ or next;
                    # Rank overlapping regions: a DIRECT region outranks a paralog one
                    # whatever their strengths, because it is evidence about this gene
                    # rather than evidence transferred onto it. Within a source, Strong
                    # outranks Moderate. (Paralog rows are capped at Moderate upstream.)
                    my $rank = (($pf[1] // "") eq "PERv1_direct" ? 2 : 0)
                             + ($st eq "Strong" ? 1 : 0);
                    next if $rank <= $pm1_rank;
                    $pm1_rank     = $rank;
                    $pm1_strength = $st;
                    $pm1_src      = $pf[1] // "";
                    # arm + region + range + stats; the strength itself is omitted
                    # because acmg_criteria already carries it.
                    $pm1_detail   = join("/", $pf[1], @pf[3 .. $#pf]);
                }
            }

            my $gt_susp = grep { /^(lowDP|lowGQ|AB_)/ } @qc;   # GT/DP suspicious?
            my ($acmg_class,$acmg_crit) = acmg_classify(
                consequence=>$consequence, lof_type=>$lof_type, loftee=>$loftee,
                freq=>$freq, nhom=>$g_nhom, revel=>$revel, am_score=>$am_score,
                eve_class=>$eve_class, cadd_num=>$cadd_num, clnsig=>$clnsig, pangolin=>$pangolin,
                ac=>$g_ac, an=>$g_an, probe=>$probe, inh=>$inheritance, gt_clean=>(!$gt_susp),
                aa_crit=>$aa_crit, aa_conflict=>$aa_conflict, clnstar=>$clnstar_n,  # PS1/PM5, PP5/BP6 star-gate
                mis_oe=>mis_oe_for($gene),                                    # PP2 (missense constraint)
                pm1=>$pm1_strength,                                            # PM1 (PERv1 region)
                de_novo_mech=>(($moi // "") =~ /\bAD\b|\bXL\b/i ? 1 : 0));   # PS2/PM6 [#6]

            # A pathogenic-leaning auto-class that carries a hard benign line (ClinVar
            # B/LB >=1 star, or a frequency/homozygote criterion) is a contradiction the
            # combining rules cannot express: "Conflicting" needs BOTH sides to reach a
            # 2-tier threshold independently, so one benign criterion never blocks a
            # pathogenic call and the row reaches the curator looking clean. This is a
            # triage tool, so the class is left alone and the tension is made visible
            # instead. BP4 is excluded on purpose — a computational prediction
            # disagreeing with PVS1/PM2 is routine, not a contradiction.
            push @qc, "clinvar_conflict"
                if $acmg_class =~ /athogenic/ && $acmg_crit =~ /(?:^|,)(BP6|BS1|BS2|BA1)(?:,|$)/;

            # Which PER earned PM1, so the curator can audit the call rather than
            # take "PM1_Strong" on trust.
            # Which region earned PM1, and from which arm, so the curator can tell
            # "this gene's own hotspot" from "inherited from its paralogs" at a glance.
            push @qc, "PM1:$pm1_detail"
                if $acmg_crit =~ /(?:^|,)PM1(?:_Strong)?(?:,|$)/ && $pm1_detail ne "";

            # PM1 and PS1/PM5 both drawing on ClinVar at the same residue is not
            # blocked -- PER regional evidence is validated independently of the
            # individual submissions -- but ClinGen SVI cautions against reusing one
            # piece of evidence twice, so the co-occurrence is surfaced for the
            # curator to rule on.
            push @qc, "PM1_with_" . lc($aa_crit)
                if $aa_crit && $acmg_crit =~ /(?:^|,)PM1(?:_Strong)?(?:,|$)/;
            $qc_flag = join(";", @qc);

            # Combined HGVS: TRANSCRIPT:c.… (p.…)  [protein accession stripped from HGVSp].
            my $hgvs = $hgvsc;
            if ($hgvsp ne "") {
                (my $p = $hgvsp) =~ s/^[^:]*://;     # drop ENSP…: prefix, keep p.…
                $hgvs = ($hgvs ne "") ? "$hgvs ($p)" : $p if $p ne "";
            }

            # Recessive carrier (g4e primary AR gene, MOI from panel) — dropped later
            # unless biallelic (hom or comp-het). Mirrors the ACMG SF AR rule [#6].
            # [#5] Recessive-carrier drop applies only to genes that are recessive
            # AND NOT also dominant. Dual-inheritance genes (MOI has AD too, e.g.
            # "AD, AR") keep their solitary hets as dominant candidates; the biallelic
            # recessive reading is still surfaced via the per-gene comp-het flag below.
            my $rec_ar = ($class eq "primary" && moi_recessive($moi) && !moi_dominant($moi)) ? 1 : 0;

            push @rows, {
                vid=>$my_id, gene=>$gene, zyg=>$zyg, mat=>$in_m, pat=>$in_f,
                class=>$class, sf_ar=>$sf_ar, rec_ar=>$rec_ar, qc=>$qc_flag, probe=>$probe,
                transcript=>$transcript, mane_select=>(field(\@r,$i{mane_select}) ne "" ? 1 : 0),
                data=>{
                    # END spans the REF allele (start + len(REF) - 1), so an indel reports
                    # its true footprint and the row round-trips into the 5-field
                    # chr-start-end-ref-alt form that resolve_variant accepts.
                    chr=>$chr, start=>$start, end=>$start + length($ref) - 1,
                    ref=>$ref, alt=>$alt,
                    gene=>$gene, strand=>$strand,
                    consequence=>$consequence, hgvs=>$hgvs,
                    revel=>$revel, eve_class=>$eve_class, eve_score=>$eve_score, cadd=>$cadd,
                    am_class=>$am_class, am_score=>$am_score, pangolin_score=>$pangolin,
                    clinvar_sig=>$clnsig, clinvar_stars=>$clnstar_n,
                    clinvar_disease=>$clndn, clinvar_aa=>$aa_detail,
                    loftee=>$loftee,
                    gnomAD_ac=>$g_ac, gnomAD_an=>$g_an, gnomAD_af=>sprintf("%.5f",$freq),
                    gnomAD_nhomalt=>$g_nhom, gnomAD_filter=>$g_filter,
                    zygosity=>$zyg, GT=>$gt, DP=>$dp, GQ=>$gq, AB=>$ab,
                    GT_SOURCE=>$gtsrc, NCALLERS=>$ncallers, CONF=>$conf,
                    inheritance=>$inheritance, kept_by=>$kept_by,
                    acmg_class=>$acmg_class, acmg_criteria=>$acmg_crit, flags=>"",
                    Association=>($assoc//""), MOI=>$moi, GDV=>($gdv//""),
                },
            };
        }
    }
    close $pfh;

    # ── (Re)write the Pangolin candidate input for this run (skipped in lookup) ──
    my $csv = "$proband.$PANEL_TAG.pangolin_input.csv";
    unless ($LOOKUP) {
        open my $cfh, ">", $csv or die "$csv: $!";
        print $cfh "CHROM,POS,REF,ALT\n";
        print $cfh "$emit{$_}\n" for sort keys %emit;
        close $cfh;
    }

    if (!$final) {
        printf "  EMIT: %d variants -> %s (%d candidate/SF + %d splice-discovery probe%s)\n",
               scalar(keys %emit), $csv,
               scalar(keys %emit) - $stat{probes}, $stat{probes},
               ($stat{probes} == 1 ? "" : "s");
        printf "  Run Pangolin on it (run_filtering.sh) to create %s, then re-run.\n", $tsv;
        next;
    }

    # ── [#3] One MANE row per variant ──
    # A variant can annotate onto >1 MANE transcript (MANE Select + MANE Plus
    # Clinical, or overlapping gene models — e.g. MUTYH's two transcripts). Collapse
    # to a single row per variant: prefer a panel-primary row, then MANE Select over
    # MANE Plus Clinical, then the most evidence arms, then a stable transcript order.
    # Skipped for --lookup consults, which report every annotation.
    unless ($LOOKUP) {
        my %best;
        for my $row (@rows) {
            my $arms = ($row->{data}{kept_by} && $row->{data}{kept_by} ne "none")
                       ? scalar(split /;/, $row->{data}{kept_by}) : 0;
            my $key = [ ($row->{class} // "") eq "primary" ? 1 : 0,
                        $row->{mane_select} ? 1 : 0,
                        $arms,
                        -length($row->{transcript} // "") ];
            # Keyed per (variant, GENE), not per variant. Keying on coordinates alone
            # collapsed across overlapping MANE gene models — MYH11+NDE1, HPDL+MUTYH,
            # COL4A1+COL4A2 all co-locate in a real exome — so a panel candidate could
            # silently delete a reportable ACMG-SF incidental at the same position, and
            # the surviving row's gene reassignment corrupted the other gene's comp-het
            # tally. The MANE Select / MANE Plus Clinical duplication this collapse
            # exists to remove is per-gene, so per-gene is the right unit.
            my $gkey = $row->{vid} . "\t" . ($row->{gene} // "");
            my $cur = $best{$gkey};
            $best{$gkey} = { row=>$row, key=>$key }
                if !$cur || _key_gt($key, $cur->{key});
        }
        # Restore a deterministic, genomically-sorted order (hash collapse loses it).
        @rows = sort { $a->{data}{chr} cmp $b->{data}{chr}
                       || $a->{data}{start} <=> $b->{data}{start}
                       || $a->{gene} cmp $b->{gene} }
                map { $_->{row} } values %best;
    }

    # ── [#6] Per-gene recessive logic over unique variants ──
    my %gene_var;   # gene -> vid -> {zyg, mat, pat}
    my %gene_rec;   # gene -> recessive-capable (MOI mentions AR/XLR)  [#7]
    for my $row (@rows) {
        $gene_var{$row->{gene}}{$row->{vid}} = { zyg=>$row->{zyg}, mat=>$row->{mat}, pat=>$row->{pat} };
        $gene_rec{$row->{gene}} //= moi_recessive($row->{data}{MOI});
    }
    # [#7] Homozygous / compound-het flags are a RECESSIVE concept — compute them
    # only for recessive-capable genes (AR/XLR, incl. dual AD/AR). A purely dominant
    # gene carrying two independent hets never gets a spurious CompHet? label.
    my %gene_flag;
    for my $g (keys %gene_var) {
        next unless $gene_rec{$g};
        my @vids = keys %{$gene_var{$g}};
        my $hom  = grep { $gene_var{$g}{$_}{zyg} eq "hom" } @vids;
        my $hem  = grep { $gene_var{$g}{$_}{zyg} eq "hem" } @vids;
        my @het  = grep { $gene_var{$g}{$_}{zyg} eq "het" } @vids;
        my $flag = "";
        if ($hom) {
            $flag = "HOM";                                  # homozygous → recessive
        } elsif ($hem) {
            $flag = "HEM";                                  # hemizygous → complete genotype
        } elsif (@het >= 2) {
            my $mat = grep {  $gene_var{$g}{$_}{mat} && !$gene_var{$g}{$_}{pat} } @het;
            my $pat = grep { !$gene_var{$g}{$_}{mat} &&  $gene_var{$g}{$_}{pat} } @het;
            $flag = ($mat && $pat) ? "CompHet(trans)" : "CompHet?";  # phaseable only in trio
        }
        $gene_flag{$g} = $flag;
    }
    # [#1,#2] Recessive-carrier handling. A solitary het in a PURE recessive gene
    # (rec_ar, or an ACMG-SF AR gene) that is not biallelic (hom / comp-het) is a
    # carrier. DEFAULT = DROP it (g4e reports no carriers: a single het cannot explain a
    # recessive disease, and carrier states are noise for clinical interpretation). This
    # does NOT affect true compound hets — a gene with >=2 gate-passing hets is biallelic
    # (CompHet flag) and every such row is kept regardless. OPT-IN --keep-ar-carriers /
    # KEEP_AR_CARRIERS=1 instead SURFACES the strong solitary carriers (the carrier-only
    # tier: kept iff strong-evidence — ClinVar P/LP >=1*, HC-LoF, or >=2 strong predictors
    # — AND not Benign/LB; flagged flags=carrier-only), for a targeted second-hit
    # hunt (a deep-intronic / CNV partner the exome may have missed). Dual AD/AR genes are
    # NOT rec_ar (see [#5]), so their solitary hets pass through as dominant candidates.
    unless ($LOOKUP) {
        @rows = grep {
            my $biallelic = (zyg_biallelic($_->{zyg}) || ($gene_flag{$_->{gene}} || "") =~ /CompHet/);
            my $solitary_carrier = ($_->{sf_ar} || $_->{rec_ar}) && !$biallelic;
            !$solitary_carrier      ? 1                       # non-carrier or biallelic: keep
              : !$KEEP_AR_CARRIERS  ? 0                       # DEFAULT: drop solitary carriers
              : (carrier_strong_evidence($_->{data}) && !is_benign_class($_->{data}));  # opt-in: carrier-only tier
        } @rows;
    }

    # The gene-level verdict decides which rows are KEPT (above); the label written to
    # each row must still describe THAT row. A gene flagged HOM because one variant is
    # homozygous must not stamp "HOM" onto an independent heterozygous variant in the
    # same gene — the row would contradict its own zygosity column and read to a curator
    # as a biallelic finding. Only a hom row is labelled HOM, and a CompHet label goes to
    # the het rows that constitute it.
    for my $row (@rows) {
        my $gf = $gene_flag{$row->{gene}} // "";
        my $biallelic = (zyg_biallelic($row->{zyg}) || $gf =~ /CompHet/);
        # $gf is populated for recessive-capable genes ONLY, so testing it here keeps
        # the existing invariant that a purely dominant gene never carries a recessive
        # flag, while the zygosity test keeps the label true of THIS row.
        my $own = ($gf eq "HOM"     && $row->{zyg} eq "hom") ? "HOM"
                : ($gf eq "HEM"     && $row->{zyg} eq "hem") ? "HEM"
                : ($gf =~ /CompHet/ && $row->{zyg} eq "het") ? $gf
                :                                              "";
        $row->{rec_label} =
            $own ne ""                                           ? $own
          : (($row->{sf_ar} || $row->{rec_ar}) && !$biallelic)   ? "carrier-only"
          :                                                        "";
        # ── Single consolidated `flags` column ──
        # The recessive verdict and the QC/artifact flags used to occupy two separate
        # columns 4 apart in a 43-column table, so a curator scanning left to right met
        # the interpretive flag and the technical one at different times. They are one
        # answer to one question — "is there anything about this row you should know
        # before reading it?" — so they are one column, recessive verdict first.
        $row->{data}{flags} = join(";", grep { defined && $_ ne "" }
                                        ($row->{rec_label}, $row->{qc}));
    }

    # ── Write output ──
    # A single-variant consult (-v / --lookup) writes ONLY the transposed,
    # human-readable view to Lookup.<tag>.<panel>.candidatos — one
    # "field <TAB> value" line per column, "." for empty — and echoes it to
    # stdout; that IS the deliverable. Normal cohort runs write the TSV
    # candidatos table + run summary as before.
    # The "Lookup." prefix is mandatory: -v already builds its VCF under that
    # name, but --lookup derives the tag from the input VCF's basename, which for
    # a cohort VCF is exactly the proband name. Without the prefix the transposed
    # consult would overwrite that proband's real <proband>.<panel>.candidatos.
    if ($LOOKUP) {
        my $rf = ($proband =~ /^Lookup\./)
               ? "$proband.$PANEL_TAG.candidatos"
               : "Lookup.$proband.$PANEL_TAG.candidatos";
        open OUT, ">", $rf or die "out: $!";
        my $nr = 0;
        for my $row (@rows) {
            $nr++;
            my $hdr = (@rows > 1) ? sprintf("──── variant %d ────\n", $nr) : "";
            print OUT $hdr; print STDOUT $hdr;
            for my $c (@COLS) {
                my $v = $row->{data}{$c};
                $v = "." if !defined $v || $v eq "";
                printf OUT    "%-16s\t%s\n", $c, $v;
                printf STDOUT "%-16s\t%s\n", $c, $v;
            }
            if (@rows > 1) { print OUT "\n"; print STDOUT "\n"; }
        }
        close OUT;
        if (@rows) { print "  -> $rf (", scalar(@rows), " variant row(s))\n"; }
        else       { print "  no reportable MANE annotation — re-run with --all-transcripts (empty $rf)\n"; }
    } else {
        # ── Write candidatos (primary + Incidental, distinguished by the GDV column) ──
        open OUT, ">$proband.$PANEL_TAG.candidatos" or die "out: $!";
        print OUT join("\t", @COLS), "\n";
        for my $row (@rows) {
            print OUT join("\t", map { $row->{data}{$_} // "" } @COLS), "\n";
            # Same row, prefixed with the sample it came from, for the batch table
            # written once after every proband is done.
            push @batch_rows, join("\t", $proband, map { $row->{data}{$_} // "" } @COLS);
        }
        close OUT;

        # ── [#9] Run summary ──
        my ($n_prim,$n_inc) = (0,0);
        my $n_probe_kept = 0;
        my (%by_arm,%by_inh,%by_flag);
        for my $row (@rows) {
            $n_probe_kept++ if $row->{probe};
            $row->{class} eq "incidental" ? $n_inc++ : $n_prim++;
            $by_arm{$_}++ for split /;/, $row->{data}{kept_by};
            $by_inh{$row->{data}{inheritance}}++;
            $by_flag{$row->{rec_label}}++ if ($row->{rec_label} // "") ne "";
        }
        print  "  -> $proband.$PANEL_TAG.candidatos\n";
        # structural-pass is reported as UNIQUE VARIANTS (keys %emit), matching what the
        # EMIT pass prints. $stat{structural} increments per CSQ annotation, so a variant
        # on two MANE transcripts counted twice and the two passes disagreed on identical input.
        printf "  variants: %d read | %d multiallelic-skipped | %d structural-pass | %d primary + %d incidental\n",
               $stat{lines}, $stat{multiallelic}, scalar(keys %emit), $n_prim, $n_inc;
        printf "  splice discovery: %d probe(s) scored, %d rescued by Pangolin >= %s\n",
               $stat{probes}, $n_probe_kept, $SPLICE_MIN if $SPLICE_PROBE && $stat{probes};
        printf "  cohort_artifact: %d variant(s) dropped (recurrent in >=%d%% of cohort & gnomAD-absent)\n",
               $stat{cohort_dropped}, int($COHORT_MAX_FRAC*100 + 0.5)
            if $cohort_on && !$KEEP_COHORT_ARTIFACTS && $stat{cohort_dropped};
        print  "  kept_by:     ", join(", ", map { "$_=$by_arm{$_}" } sort keys %by_arm), "\n" if %by_arm;
        print  "  inheritance: ", join(", ", map { "$_=$by_inh{$_}" } sort keys %by_inh), "\n" if %by_inh;
        print  "  recessive:   ", join(", ", map { "$_=$by_flag{$_}" } sort keys %by_flag), "\n" if %by_flag;
    }
}

# ── Batch-level table: every proband's candidates in one file ──
# Written only for a normal (non-consult) run that actually produced final-pass rows.
# Column 1 is the sample; the remaining columns are exactly the per-proband table, so
# the two stay diffable and anything that reads one can read the other.
if (!$LOOKUP && @batch_rows) {
    my $bf = "batch.$PANEL_TAG.candidatos";
    open my $bfh, ">", $bf or die "$bf: $!";
    print $bfh join("\t", "sample", @COLS), "\n";
    print $bfh "$_\n" for @batch_rows;
    close $bfh;
    printf "\n-> %s (%d candidate rows across %d proband(s))\n",
           $bf, scalar(@batch_rows), scalar(@probands);
}

# ── Variant-mode cleanup: drop the annotated lookup VCF, splice scratch + tmp dir ──
if (defined $LOOKUP_ANN) {
    if ($KEEP_VCF) {
        print "kept annotated VCF: $LOOKUP_ANN\n";
    } else {
        unlink $LOOKUP_ANN, "$LOOKUP_ANN.tbi", "$LOOKUP_ANN.csi", "${LOOKUP_ANN}_summary.html";
    }
    (my $lbase = $LOOKUP_ANN) =~ s/\.germline\.vep\.vcf\.gz$//;
    unlink "$lbase.$PANEL_TAG.pangolin.tsv";        # regenerable splice scratch (score is in the .candidatos output)
    system("rm", "-rf", $LOOKUP_WORK) if defined $LOOKUP_WORK && -d $LOOKUP_WORK;
}

print "\nlisto!\n";
