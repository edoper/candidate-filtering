#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

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
#  * Inclusion (rescue) gate (OR): CADD>=22, AlphaMissense path, EVE path,
#    REVEL>=0.5, Pangolin>=0.5, ClinVar P/LP, LoF (LOFTEE HC or high-impact
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
my $CADD_MIN   = 22;     # CADD PHRED rescue threshold
my $REVEL_MIN  = 0.5;    # REVEL rescue threshold (permissive; ClinGen PP3 ~0.644)
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

# Parent carrier map: chr-pos-ref-alt -> 1 if the parent carries the ALT
# (GT contains a '1'); records with 0/0 or no-call are NOT carriers.
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
        my ($gt) = parse_call($fmt,$smp);
        my $z = zygosity($gt);
        $carry{"$chr-$pos-$ref-$alt"} = 1 if $z eq "het" || $z eq "hom";
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

#############################################################################
# Discover trios / duos from filenames
#############################################################################

my @files = glob $INPUT_GLOB;
my %file_for;
for my $f (@files) { (my $s = $f) =~ s/\.germline\.vep\.vcf\.gz$//; $file_for{$s} = $f; }

my %is_parent;
for my $s (keys %file_for) {
    $is_parent{$s} = 1 if $s =~ /^(.+)([MF])$/ && exists $file_for{$1};
}

# Proband list: forced override (--proband) or filename-based auto-discovery.
my @probands;
if (@force_probands) {
    for my $p (@force_probands) {
        die "forced proband '$p' has no '$p.germline.vep.vcf.gz' in this directory\n"
            unless exists $file_for{$p};
    }
    @probands = @force_probands;
} else {
    @probands = sort grep { !$is_parent{$_} } keys %file_for;
}

print "muestras encontradas: ", join(", ", sort keys %file_for), "\n";
print "probandos: ", join(", ", @probands), (@force_probands ? " (forced override)" : ""), "\n";

# Output column order (single definition, reused for header + rows).
my @COLS = qw(
    chr start end ref alt gene strand transcript consequence hgvs.c hgvs.p tpos
    revel eve_class eve_score cadd am_class am_score pangolin_score
    clinvar_sig clinvar_stars clinvar_disease
    loftee loftee_filter loftee_flags
    gnomAD_ac gnomAD_an gnomAD_af gnomAD_nhomalt gnomAD_filter
    zygosity GT DP GQ AB
    inheritance recessive_flag kept_by Association MOI GDV
);

#############################################################################
# Process each proband
#############################################################################

foreach my $proband (@probands) {
    my $pfile = $file_for{$proband};
    my $mfile = $file_for{$proband."M"};
    my $ffile = $file_for{$proband."F"};
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

            # EMIT pass: only candidate structural-pass variants need Pangolin.
            if (!$final) {
                $emit{$my_id} = "$chr,$start,$ref,$alt"
                    if $cand_structural && !exists $emit{$my_id};
                next;
            }
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
                push @kept, "AM"       if $am_class  =~ /athogenic/;
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

            push @rows, {
                vid=>$my_id, gene=>$gene, zyg=>$zyg, mat=>$in_m, pat=>$in_f,
                class=>$class, sf_ar=>$sf_ar,
                data=>{
                    chr=>$chr, start=>$start, end=>$start, ref=>$ref, alt=>$alt,
                    gene=>$gene, strand=>$strand, transcript=>$transcript,
                    consequence=>$consequence, 'hgvs.c'=>$hgvsc, 'hgvs.p'=>$hgvsp, tpos=>$tpos,
                    revel=>$revel, eve_class=>$eve_class, eve_score=>$eve_score, cadd=>$cadd,
                    am_class=>$am_class, am_score=>$am_score, pangolin_score=>$pangolin,
                    clinvar_sig=>$clnsig, clinvar_stars=>$clnstars, clinvar_disease=>$clndn,
                    loftee=>$loftee, loftee_filter=>$lof_filt, loftee_flags=>$lof_flag,
                    gnomAD_ac=>$g_ac, gnomAD_an=>$g_an, gnomAD_af=>sprintf("%.5f",$freq),
                    gnomAD_nhomalt=>$g_nhom, gnomAD_filter=>$g_filter,
                    zygosity=>$zyg, GT=>$gt, DP=>$dp, GQ=>$gq, AB=>$ab,
                    inheritance=>$inheritance, recessive_flag=>"", kept_by=>$kept_by,
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
    # Recessive ACMG SF genes: report biallelic only (hom or comp-het); drop carriers.
    @rows = grep {
        !( $_->{sf_ar} && !($_->{zyg} eq "hom" || ($gene_flag{$_->{gene}} || "") =~ /CompHet/) )
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
