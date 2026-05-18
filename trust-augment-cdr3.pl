#!/usr/bin/env perl

# Augment _cdr3.out with per-contig coverage stats (from _final.out)
# and optional per-contig read count (from _assign.out).
#
# Usage:
#   perl trust-augment-cdr3.pl <prefix>
#     reads:  <prefix>_cdr3.out, <prefix>_final.out, [<prefix>_assign.out]
#     writes: <prefix>_cdr3.out (in place; original backed up to _cdr3.out.bak)
#
# Added columns appended to each _cdr3.out row:
#   length  min_internal_cov  mean_cov  max_cov  [read_count]
# read_count is appended only if <prefix>_assign.out exists.

use strict ;
use warnings ;

die "Usage: perl trust-augment-cdr3.pl <prefix>\n" if (@ARGV < 1) ;
my $prefix = $ARGV[0] ;

my $cdr3File   = "${prefix}_cdr3.out" ;
my $finalFile  = "${prefix}_final.out" ;
my $assignFile = "${prefix}_assign.out" ;

die "Missing $cdr3File\n"  if (! -e $cdr3File) ;
die "Missing $finalFile\n" if (! -e $finalFile) ;

# ---------- 1. Parse _final.out for per-contig coverage ----------
# Layout per contig:
#   >contig_id [extra fields]
#   <consensus seq>
#   <A counts space-separated>
#   <C counts space-separated>
#   <G counts space-separated>
#   <T counts space-separated>
my %cov ;   # contig_id => [length, min_internal_cov, mean_cov, max_cov]

open my $FF, "<", $finalFile or die "Cannot open $finalFile: $!" ;
while (my $header = <$FF>) {
    chomp $header ;
    next if ($header !~ /^>/) ;
    my $id = (split /\s+/, $header)[0] ;
    $id =~ s/^>// ;

    my $consensus = <$FF> ;
    next unless defined $consensus ;
    chomp $consensus ;
    my $len = length($consensus) ;

    my @perBase = (0) x $len ;
    my $hasPosWeight = 1 ;
    for (my $k = 0 ; $k < 4 ; ++$k) {
        my $line = <$FF> ;
        if (!defined $line) { $hasPosWeight = 0 ; last ; }
        chomp $line ;
        my @counts = split /\s+/, $line ;
        # If posWeight is missing this contig will have a different following structure;
        # in practice TRUST4 always emits 4 rows here.
        if (scalar(@counts) < $len) { $hasPosWeight = 0 ; last ; }
        for (my $j = 0 ; $j < $len ; ++$j) {
            $perBase[$j] += $counts[$j] + 0 ;
        }
    }

    if (!$hasPosWeight || $len == 0) {
        $cov{$id} = [$len, 0, 0, 0] ;
        next ;
    }

    my ($minC, $maxC, $sumC) = ($perBase[0], $perBase[0], 0) ;
    for my $c (@perBase) {
        $minC = $c if ($c < $minC) ;
        $maxC = $c if ($c > $maxC) ;
        $sumC += $c ;
    }
    my $meanC = $sumC / $len ;
    $cov{$id} = [$len, $minC, sprintf("%.2f", $meanC), $maxC] ;
}
close $FF ;

# ---------- 2. Parse _assign.out for per-contig read count (optional) ----------
my %readCnt ;
my $hasAssign = (-e $assignFile) ? 1 : 0 ;
if ($hasAssign) {
    open my $AF, "<", $assignFile or die "Cannot open $assignFile: $!" ;
    while (my $line = <$AF>) {
        chomp $line ;
        next if ($line eq "") ;
        my @f = split /\t/, $line ;
        next if (@f < 2) ;
        # Format: read_id <TAB> contig_name [<TAB> ...]
        my $contig = $f[1] ;
        $readCnt{$contig}++ ;
    }
    close $AF ;
}

# ---------- 3. Rewrite _cdr3.out with appended columns ----------
my $bak = "${cdr3File}.bak" ;
rename $cdr3File, $bak or die "Cannot back up $cdr3File: $!" ;

open my $IN,  "<", $bak     or die "Cannot read $bak: $!" ;
open my $OUT, ">", $cdr3File or die "Cannot write $cdr3File: $!" ;
while (my $line = <$IN>) {
    chomp $line ;
    next if ($line eq "") ;
    my @f = split /\t/, $line ;
    my $id = $f[0] ;

    my ($len, $minC, $meanC, $maxC) = (0, 0, 0, 0) ;
    if (exists $cov{$id}) {
        ($len, $minC, $meanC, $maxC) = @{$cov{$id}} ;
    }

    if ($hasAssign) {
        my $rc = exists $readCnt{$id} ? $readCnt{$id} : 0 ;
        print $OUT join("\t", @f, $len, $minC, $meanC, $maxC, $rc), "\n" ;
    }
    else {
        print $OUT join("\t", @f, $len, $minC, $meanC, $maxC), "\n" ;
    }
}
close $IN ;
close $OUT ;
unlink $bak ;
