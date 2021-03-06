package edu.unc.bioinf.ubu.assembly;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

public class ReAlignerOptions extends Options {
	
	private static final String INPUT_SAM = "in";
	private static final String INPUT_SAM2 = "in2";
	private static final String OUTPUT_SAM = "out";
	private static final String OUTPUT_SAM2 = "out2";
	private static final String REFERENCE = "ref";
	private static final String REFERENCE_DIR = "ref-dir";
	private static final String TARGET_REGIONS = "targets";
	private static final String WORKING_DIR = "working";
	private static final String KMER_SIZE = "kmer";
	private static final String MIN_EDGE_FREQUENCY = "mef";
	private static final String MIN_NODE_FREQUENCY = "mnf";
	private static final String MIN_EDGE_RATIO = "mer";
	private static final String MIN_CONTIG_LENGTH = "mcl";
	private static final String MAX_POTENTIAL_CONTIGS = "mpc";
	private static final String MIN_CONTIG_RATIO = "mcr";
	private static final String MIN_CONTIG_MAPQ = "mc-mapq";
	private static final String NUM_THREADS = "threads";
	private static final String MIN_UNIQUE_READS = "mur";
	private static final String ALLOWED_MISMATCHES_FROM_CONTIG = "amc";
	private static final String SKIP_UNALIGNED_ASSEMBLY = "no-unalign";
	private static final String SMALL_ALIGNER_INDEX = "small-idx";
	
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(INPUT_SAM2, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM, "Required output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM2, "Required output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Genome reference location").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE_DIR, "Directory of individually referenced chromosomes").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS, "GTF containing target regions").withRequiredArg().ofType(String.class);
            parser.accepts(WORKING_DIR, "Working directory for intermediate output").withRequiredArg().ofType(String.class);
            parser.accepts(KMER_SIZE, "Assembly kmer size").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_EDGE_FREQUENCY, "Assembly minimum edge frequency").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_NODE_FREQUENCY, "Assembly minimum node frequency").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_EDGE_RATIO, "Assembly minimum edge ratio").withRequiredArg().ofType(Double.class);
            parser.accepts(MIN_CONTIG_LENGTH, "Assembly minimum contig length").withRequiredArg().ofType(Integer.class);
            parser.accepts(MAX_POTENTIAL_CONTIGS, "Maximum number of potential contigs for a region").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_CONTIG_RATIO, "Minimum contig length as percentage of observed region length").withRequiredArg().ofType(Double.class);
            parser.accepts(NUM_THREADS, "Number of threads").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_CONTIG_MAPQ, "Minimum contig mapping quality").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_UNIQUE_READS, "Minimum number of unique reads per node").withRequiredArg().ofType(Integer.class);
            parser.accepts(ALLOWED_MISMATCHES_FROM_CONTIG, "Maximum number of allowed mismatches from contig for a read (default zero).").withRequiredArg().ofType(Integer.class);
            parser.accepts(SKIP_UNALIGNED_ASSEMBLY, "Skip assembly of reads that do not align to contigs.");
            parser.accepts(SMALL_ALIGNER_INDEX, "Use small aligner index");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
		isValid = true;
		
		if (!getOptions().hasArgument(INPUT_SAM)) {
			isValid = false;
			System.out.println("Missing required input SAM/BAM file");
		}

		if (!getOptions().hasArgument(OUTPUT_SAM)) {
			isValid = false;
			System.out.println("Missing required input SAM/BAM file");
		}
		
		if ((!getOptions().hasArgument(REFERENCE)) && (!getOptions().hasArgument(REFERENCE_DIR))) {
			isValid = false;
			System.out.println("Missing required reference");
		}

		if ((getOptions().hasArgument(REFERENCE)) && (getOptions().hasArgument(REFERENCE_DIR))) {
			isValid = false;
			System.out.println("Please only specify one of: " + REFERENCE + ", " + REFERENCE_DIR);
		}
		
		if (!getOptions().hasArgument(TARGET_REGIONS)) {
			isValid = false;
			System.out.println("Missing required target regions");
		}
		
		if (!getOptions().hasArgument(WORKING_DIR)) {
			isValid = false;
			System.out.println("Missing required working directory");
		}
		
		if ((Integer) getOptions().valueOf(NUM_THREADS) < 1) {
			isValid = false;
			System.out.println("Num threads must be greater than zero.");
		}
		
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT_SAM);
	}
	
	public String getOutputFile() {
		return (String) getOptions().valueOf(OUTPUT_SAM);
	}
	
	public String getInputFile2() {
		return (String) getOptions().valueOf(INPUT_SAM2);
	}
	
	public String getOutputFile2() {
		return (String) getOptions().valueOf(OUTPUT_SAM2);
	}
	
	public String getReference() {
		return (String) getOptions().valueOf(REFERENCE);
	}
	
	public String getReferenceDir() {
		return (String) getOptions().valueOf(REFERENCE_DIR);
	}
	
	public String getTargetRegionFile() {
		return (String) getOptions().valueOf(TARGET_REGIONS);
	}
	
	public String getWorkingDir() {
		return (String) getOptions().valueOf(WORKING_DIR);
	}
	
	public int getKmerSize() {
		return (Integer) getOptions().valueOf(KMER_SIZE);
	}
	
	public int getMinEdgeFrequency() {
		return (Integer) getOptions().valueOf(MIN_EDGE_FREQUENCY);
	}
	
	public int getMinNodeFrequency() {
		return (Integer) getOptions().valueOf(MIN_NODE_FREQUENCY);
	}
	
	public double getMinEdgeRatio() {
		return (Double) getOptions().valueOf(MIN_EDGE_RATIO);
	}
	
	public int getMinContigLength() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_LENGTH);
	}
	
	public int getMaxPotentialContigs() {
		return (Integer) getOptions().valueOf(MAX_POTENTIAL_CONTIGS);
	}
	
	public double getMinContigRatio() {
		return (Double) getOptions().valueOf(MIN_CONTIG_RATIO);
	}
	
	public int getNumThreads() {
		return (Integer) getOptions().valueOf(NUM_THREADS);
	}
	
	public int getMinContigMapq() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_MAPQ);
	}
	
	public int getMinUniqueReads() {
		return (Integer) getOptions().valueOf(MIN_UNIQUE_READS);
	}
	
	public boolean isSkipUnalignedAssembly() {
		return getOptions().has(SKIP_UNALIGNED_ASSEMBLY);
	}
	
	public boolean useSmallAlignerIndex() {
		return getOptions().has(SMALL_ALIGNER_INDEX);
	}
	
	public boolean isValid() {
		return isValid;
	}
	
	public int getAllowedMismatchesFromContig() {
		return (Integer) getOptions().valueOf(ALLOWED_MISMATCHES_FROM_CONTIG);
	}
}
