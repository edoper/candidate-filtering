#!/usr/bin/env perl
use strict;
use warnings;

# Parse a Pangolin CSV output into a score map for filtering_r.pl.
#
#   Usage: parse_pangolin.pl <proband.pangolin.csv>  >  <proband>.pangolin.tsv
#
# Pangolin appends a 'Pangolin' column of the form
#   gene|pos:score_change|pos:score_change|warnings,gene2|...
# whose value itself contains commas. The input columns are exactly
# CHROM,POS,REF,ALT, so we split off the first 4 comma fields and treat the
# remainder as the score blob, then take the maximum |score_change| over all
# 'pos:score' tokens. Output: "chr-pos-ref-alt <tab> max_abs_delta".

my $in = shift or die "usage: $0 <pangolin_output.csv>\n";
open my $fh, "<", $in or die "$in: $!";
my $header = <$fh>;   # discard header line
while (my $l = <$fh>) {
    chomp $l;
    my @f = split /,/, $l, 5;
    next if @f < 4;
    my ($chr, $pos, $ref, $alt) = @f[0 .. 3];
    my $blob = defined $f[4] ? $f[4] : "";
    # Capture the number after each colon (the score change), skip Warnings:text.
    my @nums = $blob =~ /:\s*(-?\d+(?:\.\d+)?)/g;
    next unless @nums;                      # variant not in a gene / no score
    my $max = 0;
    for my $n (@nums) { my $a = abs($n); $max = $a if $a > $max; }
    printf "%s-%s-%s-%s\t%.2f\n", $chr, $pos, $ref, $alt, $max;
}
close $fh;
