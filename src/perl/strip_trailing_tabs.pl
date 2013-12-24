#!/usr/bin/perl

use strict;
use Getopt::Long;

# This script exists to workaround the quote issue in Seqware's FTL's

my $input;
my $temp;

GetOptions("input=s"=>\$input, "temp-file=s"=>\$temp);

system("mv $input $temp") == 0 or die $!;
my $command = "sed 's/\\t\$//g' $temp > $input";
print "$command\n";
system($command) == 0 or die $!;