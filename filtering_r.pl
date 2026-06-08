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
#  * Reads *.germline.vep.vcf.gz directly; auto-discovers families by filename:
#    <base>=proband, <base>M=mother, <base>F=father.
#  * Multiallelic input is assumed pre-split (vep_annotate.sh runs
#    `bcftools norm -m-any`). Any residual multiallelic row is skipped+counted.
#                                                                      [#2]
#  * Structural gates (AND): MANE transcript + whitelisted consequence
#    (matched per-&-atom) + panel gene + rare. Consequence is split on '&' so
#    unanticipated compound terms are not silently dropped.            [#1]
#  * Frequency threshold is MOI-aware: recessive genes tolerate a higher
#    gnomAD AF than dominant genes.                                    [#6]
#  * Inclusion (rescue) gate (OR): CADD>=25.3, AlphaMissense>=0.792, EVE path,
#    REVEL>=0.644, Pangolin>=0.5, ClinVar P/LP, LoF (LOFTEE HC or high-impact
#    truncating consequence). Each surviving row records which arm(s) fired in
#    a `kept_by` column.                                          [#4,#8]
#  * Genotype-aware: proband zygosity / DP / GQ / allele-balance columns from
#    FORMAT; parent inheritance uses parental GT (carrier = non-ref), not mere
#    site presence.                                                [#5,#10]
#  * Per-gene recessive logic: homozygous and compound-het (trans where
#    phaseable) flags.                                                 [#6]
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
#    supporting-benign, as 2015 has no benign-Moderate). Other criteria: PVS1,
#    PS2/PM6, PM2/PM4, PP5, BA1/BS1/BS2/BP6/BP7.                            [#2]
#  * QC / artifact flags (qc_flag): lowDP, lowGQ, AB_het/AB_hom, homopolymer
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
#############################################################################

# ── Tunable thresholds (single source of truth) ──
my $INPUT_GLOB = '*.germline.vep.vcf.gz';
my $FREQ_AD    = 0.01;   # max gnomAD AF (%) for dominant genes  (1 in 10,000)
my $FREQ_AR    = 1.0;    # max gnomAD AF (%) for recessive genes (carrier freq)
my $CADD_MIN   = 25.3;   # CADD PHRED rescue threshold
my $REVEL_MIN  = 0.644;  # REVEL rescue threshold (ClinGen PP3)
my $AM_MIN     = 0.792;  # AlphaMissense pathogenicity rescue threshold
my $SPLICE_MIN = 0.5;    # Pangolin |delta| splice rescue threshold

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
my $REF_FASTA   = '/home/edo/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa';
my $HAVE_REF    = -e "$REF_FASTA.fai";   # samtools-indexed reference for homopolymer check

# ── Automated ACMG/AMP classification (InterVar-style, triage only) [#2] ──
my $PM2_AC_MAX  = 1;       # gnomAD AC at/below -> PM2 (absent=0 or singleton=1)
my $BS1_FREQ    = 1.0;     # gnomAD AF (%) at/above -> BS1 (too common for rare disease)
my $BA1_FREQ    = 5.0;     # gnomAD AF (%) at/above -> BA1 (benign standalone)
my $BS2_NHOM    = 10;      # gnomAD homozygotes at/above -> BS2
my $BP4_REVEL   = 0.290;   # REVEL at/below -> BP4 (computational benign, ClinGen)

#############################################################################
# Reference hashes
#############################################################################

open MANE, "<mane-plus-clinical-names.txt" or die "mane: $!";
my %mane;
while (my $m = <MANE>) { chomp $m; my ($a,$b) = split /\t/, $m; $mane{$a} = $b; }
close MANE;
print "hash mane, listo!\n";

# ── Argument parsing ──
#   --proband NAME   (repeatable) force NAME as a proband, overriding filename-
#                    based auto-discovery. NAME must have <NAME>.germline.vep.vcf.gz.
#   <positional>     genes-of-interest file (gene panel override; see below).
#
# Gene panel: gene -> "Association \t MOI \t GDV".
#   Default source = g4e-2025.txt (4 columns). A genes-of-interest file (one gene
#   symbol per line; plain symbols -> Association/MOI/GDV = "NA", or full 4-column
#   g4e format) overrides it. '#' comments and blanks are skipped.
my (@force_probands, $GENES_FILE);
while (@ARGV) {
    my $a = shift @ARGV;
    if    ($a eq '--proband' || $a eq '-p') { push @force_probands, (shift @ARGV // ''); }
    elsif ($a =~ /^--proband=(.+)$/)        { push @force_probands, $1; }
    else                                    { $GENES_FILE = $a; }   # positional = genes file
}
my $PANEL = (defined $GENES_FILE && $GENES_FILE ne "") ? $GENES_FILE : "g4e-2025.txt";
my $custom_panel = (defined $GENES_FILE && $GENES_FILE ne "") ? 1 : 0;

# Output tag = panel basename without extension (e.g. g4e-2025, Hyperparathyroidism).
# All per-run outputs are namespaced by it so different panels don't overwrite.
my $PANEL_TAG = $PANEL;
$PANEL_TAG =~ s{.*/}{};
$PANEL_TAG =~ s/\.[^.]+$//;

open PANEL, "<", $PANEL or die "gene panel '$PANEL': $!";
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

# Consequence whitelist (atomic terms recommended; compound entries harmless).
open VAR, "<typevar.txt" or die "typevar: $!";
my %varfilter;
while (my $t = <VAR>) { chomp $t; my ($c,$d) = split /\t/, $t; $varfilter{$c} = $d//""; }
close VAR;
print "hash var, listo!\n";

# ACMG SF v3.2 genes: gene -> "condition \t MOI \t report_category". Always loaded
# (independent of the candidate panel). Non-fatal if absent.
my %acmg;
if (open my $afh, "<", $ACMG_FILE) {
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

# Zygosity from a GT string: hom (alt/alt), het (ref/alt), ref, or "" (no-call).
sub zygosity {
    my ($gt) = @_;
    return "" unless defined $gt && $gt ne "";
    my @a = split /[\/|]/, $gt;
    return "" unless @a >= 2 && $a[0] ne "." && $a[1] ne ".";
    my $n1 = grep { $_ eq "1" } @a;
    return $n1 >= 2 ? "hom" : $n1 == 1 ? "het" : "ref";
}

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
        $carry{"$chr-$pos-$ref-$alt"} = "$gt:$dp:$gq" if $z eq "het" || $z eq "hom";
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
# ClinVar review status (CLNREVSTAT) -> star count (0-4).
sub clinvar_stars {
    my ($s) = @_;
    return 0 unless defined $s && $s ne "";
    return 4 if $s =~ /practice_guideline/i;
    return 3 if $s =~ /expert_panel/i;
    return 2 if $s =~ /multiple_submitters/i;
    return 1 if $s =~ /single_submitter|criteria_provided/i;
    return 0;
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
    push @P, "PVS1" if $v{loftee} eq "HC" || ($v{lof_type} && $v{loftee} ne "LC");
    if ($v{inh} eq "DN") {                        # trio de novo (relatedness assumed)
        push @P, ($v{gt_clean} ? "PS2" : "PM6");
    } elsif ($v{inh} =~ m{^DN/} && $v{de_novo_mech}) { push @P, "PM6"; }
    push @P, "PM4" if $v{consequence} =~ /inframe_(insertion|deletion)|stop_lost/;
    push @P, "PM2" if $v{ac} ne "" && $v{ac} <= $PM2_AC_MAX;       # absent or singleton
    push @P, "PP5" if clinvar_pathogenic($v{clnsig});

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

    # Benign criteria
    push @B, "BA1" if $v{freq} >= $BA1_FREQ;
    push @B, "BS1" if $v{freq} >= $BS1_FREQ && $v{freq} < $BA1_FREQ;
    push @B, "BS2" if $v{nhom} ne "" && $v{nhom} >= $BS2_NHOM;
    push @B, "BP6" if clinvar_benign($v{clnsig});
    push @B, "BP7" if $v{consequence} =~ /synonymous_variant/
                   && (($v{pangolin} eq "" ? 0 : $v{pangolin}) < 0.2);

    # ── Combine per ACMG 2015. Count by strength tier; graded PP3/BP4 contribute
    #    at their tier (PP3_Strong->PS, _Moderate->PM, _Supporting->PP;
    #    BP4_Strong->BS, BP4_Moderate/_Supporting->BP). ──
    my $pvs = grep { $_ eq "PVS1" } @P;
    my $ps  = grep { $_ eq "PS2"  } @P;
    my $pm  = grep { /^PM\d/      } @P;            # PM2, PM4, PM6
    my $pp  = grep { $_ eq "PP5"  } @P;            # PP5
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
}

print "muestras encontradas: ", join(", ", sort keys %file_for), "\n";
print "probandos: ", join(", ", @probands), (@force_probands ? " (forced override)" : ""), "\n";

# Output column order (single definition, reused for header + rows).
my @COLS = qw(
    chr start end ref alt gene strand consequence hgvs
    revel eve_class eve_score cadd am_class am_score pangolin_score
    clinvar_sig clinvar_disease
    loftee
    gnomAD_ac gnomAD_an gnomAD_af gnomAD_nhomalt gnomAD_filter
    zygosity GT DP GQ AB
    inheritance recessive_flag kept_by acmg_class acmg_criteria qc_flag
    Association MOI GDV
);

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
    my $final = -e $tsv;
    my $pscore = $final ? load_scores($tsv) : {};

    my $kind = ($have_m && $have_f) ? "trio" : ($have_m||$have_f) ? "duo" : "singleton";
    print "\ntrabajando con $proband ($kind, ",
          ($final ? "FINAL: scores from $tsv" : "EMIT: no scores yet"), ")";
    print " | madre: $mfile" if $have_m;
    print " | padre: $ffile" if $have_f;
    print "\n";

    # Resolve CSQ indices for this proband's layout.
    my $col = csq_columns($pfile);
    my %i = (
        gene          => resolve($col,'SYMBOL'),
        strand        => resolve($col,'STRAND'),
        transcript    => resolve($col,'Feature'),
        consequence   => resolve($col,'Consequence'),
        hgvsc         => resolve($col,'HGVSc'),
        hgvsp         => resolve($col,'HGVSp'),
        tpos          => resolve($col,'cDNA_position'),
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
    my %stat = (lines=>0, multiallelic=>0, structural=>0);

    my %emit;            # EMIT pass: unique candidate variants
    my @rows;            # FINAL pass: buffered rows (for per-gene recessive logic)

    my $pfh = open_vcf($pfile);
    while (my $v = <$pfh>) {
        chomp $v;
        next if $v =~ /^#/;
        $stat{lines}++;

        my @c = split /\t/, $v;
        my ($chr,$start,$ref,$alt,$info,$fmt,$smp) = @c[0,1,3,4,7,8,9];

        if ($alt =~ /,/) { $stat{multiallelic}++; next; }   # [#2] should be pre-split

        my ($csq) = $info =~ /(?:^|;)CSQ=([^;]*)/;
        next unless defined $csq;

        my $my_id = "$chr-$start-$ref-$alt";

        # Proband genotype (same for all transcripts of this variant) [#5].
        my ($gt,$dp,$gq,$adr,$ada) = parse_call($fmt,$smp);
        my $zyg = zygosity($gt);
        my $ab  = ($ada ne "" && ($adr+$ada) > 0) ? sprintf("%.2f", $ada/($adr+$ada)) : "";

        foreach my $fila (split /,/, $csq) {
            my @r = split /\|/, $fila, -1;
            my $gene        = field(\@r,$i{gene});
            my $transcript  = field(\@r,$i{transcript});
            my $consequence = field(\@r,$i{consequence});

            next unless exists $mane{$transcript};   # both paths require MANE
            my $in_panel = exists $epigenes{$gene};
            my $in_acmg  = exists $acmg{$gene};
            next unless $in_panel || $in_acmg;

            my $g_ac = field(\@r,$i{g_ac}); $g_ac = ($g_ac eq "") ? 0 : $g_ac;
            my $g_an = field(\@r,$i{g_an}); $g_an = ($g_an eq "") ? 0 : $g_an;
            my $freq = ($g_an > 0) ? ($g_ac/$g_an)*100 : 0;

            # Candidate structural gate: panel gene + whitelisted consequence + rare
            # (MOI-aware: recessive genes tolerate higher carrier freq). [#1,#6]
            my $cand_structural = 0;
            my ($p_assoc,$p_moi,$p_gdv) = ("","","");
            if ($in_panel && consequence_ok($consequence)) {
                ($p_assoc,$p_moi,$p_gdv) = split /\t/, $epigenes{$gene};
                $p_moi //= "";
                my $recessive = ($p_moi =~ /\bAR\b|XLR|recessiv/i) ? 1 : 0;
                $cand_structural = ($freq <= ($recessive ? $FREQ_AR : $FREQ_AD)) ? 1 : 0;
            }

            # Collect candidate structural-pass variants for Pangolin. Done in BOTH
            # passes so the input CSV is (re)written with the full set even in the
            # final pass (otherwise it would be emptied and the score cache wiped).
            $emit{$my_id} = "$chr,$start,$ref,$alt"
                if $cand_structural && !exists $emit{$my_id};
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
            my $clnsig    = field(\@r,$i{clnsig});
            my $clnstars  = field(\@r,$i{clnstars});
            my $clndn     = field(\@r,$i{clndn});
            my $strand    = field(\@r,$i{strand});
            my $hgvsc     = field(\@r,$i{hgvsc});
            my $hgvsp     = field(\@r,$i{hgvsp});
            my $tpos      = field(\@r,$i{tpos});
            my $pangolin  = exists $pscore->{$my_id} ? $pscore->{$my_id} : "";
            my $cadd_num  = ($cadd eq "") ? 0 : $cadd;
            my $lof_type  = grep { $LOF_CONS{$_} } split /&/, $consequence;

            # ── Primary candidate inclusion gate (OR) [#4,#8] ──
            my (@kept, $class, $assoc, $moi, $gdv);
            if ($cand_structural) {
                push @kept, "CADD"     if $cadd_num >= $CADD_MIN;
                push @kept, "AM"       if $am_score ne "" && $am_score >= $AM_MIN;
                push @kept, "EVE"      if $eve_class =~ /athogenic/;
                push @kept, "REVEL"    if $revel ne "" && $revel >= $REVEL_MIN;
                push @kept, "Pangolin" if $pangolin ne "" && $pangolin >= $SPLICE_MIN;
                push @kept, "ClinVar"  if clinvar_pathogenic($clnsig);
                push @kept, "LoF"      if $loftee eq "HC" || ($lof_type && $loftee ne "LC");
                if (@kept) { $class = "primary"; ($assoc,$moi,$gdv) = ($p_assoc,$p_moi,$p_gdv); }
            }

            # ── Incidental (ACMG SF) — stringent; only if not a primary candidate ──
            # ClinVar P/LP (>=1 star) reported regardless of frequency (known founder
            # alleles); novel LoF / >=2-strong-computational tiers require rarity.
            my $sf_ar = 0;
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

            next unless @kept;
            my $kept_by = join(";", @kept);

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
                push @qc, "AB_hom" if $zyg eq "hom" && $ab < 0.85;
            }
            push @qc, "homopolymer"   if homopolymer_context($chr,$start,$ref,$alt);
            push @qc, "inh_lowqual"   if $inh_lowqual;
            push @qc, "DN_unconfirmed" if $inheritance =~ /^DN/;   # no parental ref depth
            my $qc_flag = join(";", @qc);

            # ── [#2] Automated ACMG/AMP classification (triage) ──
            my $gt_susp = grep { /^(lowDP|lowGQ|AB_)/ } @qc;   # GT/DP suspicious?
            my ($acmg_class,$acmg_crit) = acmg_classify(
                consequence=>$consequence, lof_type=>$lof_type, loftee=>$loftee,
                freq=>$freq, nhom=>$g_nhom, revel=>$revel, am_score=>$am_score,
                eve_class=>$eve_class, cadd_num=>$cadd_num, clnsig=>$clnsig, pangolin=>$pangolin,
                ac=>$g_ac, inh=>$inheritance, gt_clean=>(!$gt_susp),
                de_novo_mech=>(($moi // "") =~ /\bAD\b|\bXL\b/i ? 1 : 0));   # PS2/PM6 [#6]

            # Combined HGVS: TRANSCRIPT:c.… (p.…)  [protein accession stripped from HGVSp].
            my $hgvs = $hgvsc;
            if ($hgvsp ne "") {
                (my $p = $hgvsp) =~ s/^[^:]*://;     # drop ENSP…: prefix, keep p.…
                $hgvs = ($hgvs ne "") ? "$hgvs ($p)" : $p if $p ne "";
            }

            # Recessive carrier (g4e primary AR gene, MOI from panel) — dropped later
            # unless biallelic (hom or comp-het). Mirrors the ACMG SF AR rule [#6].
            my $rec_ar = ($class eq "primary" && (($moi // "") =~ /\bAR\b|XLR|recessiv/i)) ? 1 : 0;

            push @rows, {
                vid=>$my_id, gene=>$gene, zyg=>$zyg, mat=>$in_m, pat=>$in_f,
                class=>$class, sf_ar=>$sf_ar, rec_ar=>$rec_ar,
                data=>{
                    chr=>$chr, start=>$start, end=>$start, ref=>$ref, alt=>$alt,
                    gene=>$gene, strand=>$strand,
                    consequence=>$consequence, hgvs=>$hgvs,
                    revel=>$revel, eve_class=>$eve_class, eve_score=>$eve_score, cadd=>$cadd,
                    am_class=>$am_class, am_score=>$am_score, pangolin_score=>$pangolin,
                    clinvar_sig=>$clnsig, clinvar_disease=>$clndn,
                    loftee=>$loftee,
                    gnomAD_ac=>$g_ac, gnomAD_an=>$g_an, gnomAD_af=>sprintf("%.5f",$freq),
                    gnomAD_nhomalt=>$g_nhom, gnomAD_filter=>$g_filter,
                    zygosity=>$zyg, GT=>$gt, DP=>$dp, GQ=>$gq, AB=>$ab,
                    inheritance=>$inheritance, recessive_flag=>"", kept_by=>$kept_by,
                    acmg_class=>$acmg_class, acmg_criteria=>$acmg_crit, qc_flag=>$qc_flag,
                    Association=>($assoc//""), MOI=>$moi, GDV=>($gdv//""),
                },
            };
        }
    }
    close $pfh;

    # ── Always (re)write the Pangolin candidate input for this run ──
    my $csv = "$proband.$PANEL_TAG.pangolin_input.csv";
    open my $cfh, ">", $csv or die "$csv: $!";
    print $cfh "CHROM,POS,REF,ALT\n";
    print $cfh "$emit{$_}\n" for sort keys %emit;
    close $cfh;

    if (!$final) {
        printf "  EMIT: %d structural-pass variants -> %s\n", scalar(keys %emit), $csv;
        printf "  Run Pangolin on it (run_filtering.sh) to create %s, then re-run.\n", $tsv;
        next;
    }

    # ── [#6] Per-gene recessive logic over unique variants ──
    my %gene_var;   # gene -> vid -> {zyg, mat, pat}
    for my $row (@rows) {
        $gene_var{$row->{gene}}{$row->{vid}} = { zyg=>$row->{zyg}, mat=>$row->{mat}, pat=>$row->{pat} };
    }
    my %gene_flag;
    for my $g (keys %gene_var) {
        my @vids = keys %{$gene_var{$g}};
        my $hom  = grep { $gene_var{$g}{$_}{zyg} eq "hom" } @vids;
        my @het  = grep { $gene_var{$g}{$_}{zyg} eq "het" } @vids;
        my $flag = "";
        if ($hom) {
            $flag = "HOM";                                  # homozygous → recessive
        } elsif (@het >= 2) {
            my $mat = grep {  $gene_var{$g}{$_}{mat} && !$gene_var{$g}{$_}{pat} } @het;
            my $pat = grep { !$gene_var{$g}{$_}{mat} &&  $gene_var{$g}{$_}{pat} } @het;
            $flag = ($mat && $pat) ? "CompHet(trans)" : "CompHet?";  # phaseable only in trio
        }
        $gene_flag{$g} = $flag;
    }
    # Recessive genes (g4e primary AR + ACMG SF AR): report biallelic only (hom or
    # comp-het); drop solitary carriers. g4e wants no carriers at all [#6].
    @rows = grep {
        my $biallelic = ($_->{zyg} eq "hom" || ($gene_flag{$_->{gene}} || "") =~ /CompHet/);
        !( ($_->{sf_ar} || $_->{rec_ar}) && !$biallelic )
    } @rows;

    for my $row (@rows) {
        $row->{data}{recessive_flag} =
            ($row->{zyg} eq "hom") ? "HOM" : ($gene_flag{$row->{gene}} || "");
    }

    # ── Write candidatos (primary + Incidental, distinguished by the GDV column) ──
    open OUT, ">$proband.$PANEL_TAG.candidatos" or die "out: $!";
    print OUT join("\t", @COLS), "\n";
    for my $row (@rows) {
        print OUT join("\t", map { $row->{data}{$_} // "" } @COLS), "\n";
    }
    close OUT;

    # ── [#9] Run summary ──
    my ($n_prim,$n_inc) = (0,0);
    my (%by_arm,%by_inh,%by_flag);
    for my $row (@rows) {
        $row->{class} eq "incidental" ? $n_inc++ : $n_prim++;
        $by_arm{$_}++ for split /;/, $row->{data}{kept_by};
        $by_inh{$row->{data}{inheritance}}++;
        $by_flag{$row->{data}{recessive_flag}}++ if $row->{data}{recessive_flag} ne "";
    }
    print  "  -> $proband.$PANEL_TAG.candidatos\n";
    printf "  variants: %d read | %d multiallelic-skipped | %d structural-pass | %d primary + %d incidental\n",
           $stat{lines}, $stat{multiallelic}, $stat{structural}, $n_prim, $n_inc;
    print  "  kept_by:     ", join(", ", map { "$_=$by_arm{$_}" } sort keys %by_arm), "\n" if %by_arm;
    print  "  inheritance: ", join(", ", map { "$_=$by_inh{$_}" } sort keys %by_inh), "\n" if %by_inh;
    print  "  recessive:   ", join(", ", map { "$_=$by_flag{$_}" } sort keys %by_flag), "\n" if %by_flag;
}

print "\nlisto!\n";
