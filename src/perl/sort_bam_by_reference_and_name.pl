#
# Sorts a bam file first by reference then by name
#
use strict;

use File::Path;
use Getopt::Long;
use Sys::Hostname;

my $getOptResult = GetOptions(
	'input=s'    => \( my $input_bam     ),
	'output=s'   => \( my $output_bam    ),
	'temp-dir=s' => \( my $temp_dir      ),
	'samtools=s' => \( my $samtools_exec ),	
);

#my $input_bam = "/home/lisle/gtof/conv/really_small_sorted_by_coord.bam";
#my $output_bam = "/home/lisle/gtof/conv/really_small_sorted_by_ref_and_read.bam";

if (!(defined($input_bam))) {
	usage();
} 

print "Starting " . scalar(localtime) . "\n";

#Create Temp dir: <pid>_<hostname>_<time>
#my $tdir = "/home/lisle/gtof/" . "$$" . "_" . hostname . "_" . time;

my $tdir = $temp_dir . "$$" . "_" . hostname . "_" . time;

print "Using temp dir: $tdir\n";

mkpath($tdir, { mode => 0775 }) or die "$! Unable to create $tdir"; 

# Get sequence references from sam file
my $references_str = `samtools view -H $input_bam | grep \@SQ | cut -f 2 | cut -d : -f 2` or die $!;

my @references = split(/\n/,$references_str);

my $sorted_files = "";

# Loop through each reference creating sorted by name bam specific to the reference
foreach (@references) {
	my $reference = $_;
	print scalar(localtime) . " Sorting $reference\n";
	my $sorted_file = "$tdir/$reference";
	my $command = "samtools view -b $input_bam $reference | samtools sort -n - $sorted_file";
	print "$command\n"; 
	my $ret = system($command);
	print "retVal: $ret\n"; 
	
	$sorted_files = $sorted_files . " " . $sorted_file . ".bam";
}

# Concatenate the sorted bam files
my $command = "samtools cat -o $output_bam $sorted_files";
print "$command\n";
my $ret = system($command);
print "last retVal: $ret\n";

print "Done " . scalar(localtime) . "\n";

sub usage() {
	print "perl sort_bam_by_reference_and_name.pl --input <sorted and indexed input bam file> --output <output bam file> --temp-dir <temp directory> --samtools <path to samtools>\n";
	exit(1);
}