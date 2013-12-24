#!/usr/bin/perl

use strict;
use Getopt::Long;

my $out = '-';
my $q = 75;
my @col;
my @also;
my $names = 1;
my $target = 1000;
my $skip = 0;
my $min=1;
GetOptions("quant=i"=>\$q, "target=i"=>\$target, "col=i@"=>\@col, "out=s"=>\$out, "also=i@"=>\@also, "skip=i"=>\$skip, "min=i"=>\$min);

my $in = shift @ARGV;

die usage() unless $in && @col;

open(OUT, ($out eq '-') ? '<&STDOUT' : ">$out") || die "Can't open $out\n";
open(IN, ($in eq '-') ? '<&STDIN' : $in) || die "Can't open $in\n";

@also = (1) if !@also && !grep {$_ eq '1'} @col;

map {$_--} @col;
map {$_--} @also;

my @d;
my $cnt = 0;
my $head ='';
while(<IN>) {
        if ($skip) {
                --$skip;
                $head .= $_;
                next;
        }
        chomp;
        my @f = split /\t/;
        if ($col[0] eq '-2') {
                @col = (1..$#f);
        }
        for (@col) {
                push @{$d[$_]}, $f[$_];
        }
        for (@also) {
                push @{$d[$_]}, $f[$_];
        }
        ++$cnt;
}
for (@col) {
        my @t = grep {$_>=$min} @{$d[$_]};
        @t = sort {$a <=> $b} @t;
        my $t=quantile(\@t, $q/100);
        for (@{$d[$_]}) {
                $_= sprintf "%.4f", $target*$_/$t;
        }
}

my @out = (sort {$a <=> $b} (@col, @also));

print OUT $head;

for (my $i=0;$i<$cnt;++$i) {
        for my $j (@out) {
                print OUT "\t" unless $j == $out[0];
                print OUT $d[$j][$i];
        }
        print OUT "\n";
}


sub usage {
<<EOF;
Usage: $0 -c COL [opts] FILE

Returns an upper quartile normalization of data in column(s) COL
of file FILE.

Col is 1-based, zeroes are ignores when calculating upper quartile

Options:
   -c|col COL    normalize this column of data (can specify more than once, or -1 for all but first col)
   -q|quant INT  quantile to use (75)
   -t|target INT target to use (1000)
   -a|also COL   output these columns also
   -o|out FILE   output to this file instead of stdout
   -m|min INT    minimum value (1)
   -s|skip INT   skip header rows
EOF
}

sub quantile {
        my ($a,$p) = @_;
        my $l = scalar(@{$a});
        my $t = ($l-1)*$p;
        my $v=$a->[int($t)];
        if ($t > int($t)) {
                return $v + $p * ($a->[int($t)+1] - $v);
        } else {
                return $v;
        }
}
