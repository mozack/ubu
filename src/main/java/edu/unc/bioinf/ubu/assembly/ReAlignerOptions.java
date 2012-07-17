package edu.unc.bioinf.ubu.assembly;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

public class ReAlignerOptions extends Options {
	
	private static final String INPUT_SAM = "in";
	private static final String OUTPUT_SAM = "out";
	private static final String REFERENCE = "ref";
	private static final String TARGET_REGIONS = "targets";
	private static final String WORKING_DIR = "working";
	private static final String KMER_LENGTH = "kmer";
	private static final String MIN_EDGE_FREQUENCY = "mef";
	private static final String MIN_NODE_FREQUENCY = "mnf";
	private static final String MIN_EDGE_RATIO = "mer";
	private static final String MIN_CONTIG_LENGTH = "mcl";
	
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM, "Required output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Genome reference location").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS, "GTF containing target regions").withRequiredArg().ofType(String.class);
            parser.accepts(WORKING_DIR, "Working directory for intermediate output").withRequiredArg().ofType(String.class);            
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
		
		if (!getOptions().hasArgument(REFERENCE)) {
			isValid = false;
			System.out.println("Missing required reference");
		}
		
		if (!getOptions().hasArgument(TARGET_REGIONS)) {
			isValid = false;
			System.out.println("Missing required target regions");
		}
		
		if (!getOptions().hasArgument(WORKING_DIR)) {
			isValid = false;
			System.out.println("Missing required working directory");
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
	
	public String getReference() {
		return (String) getOptions().valueOf(REFERENCE);
	}
	
	public String getTargetRegionFile() {
		return (String) getOptions().valueOf(TARGET_REGIONS);
	}
	
	public String getWorkingDir() {
		return (String) getOptions().valueOf(WORKING_DIR);
	}
	
	public boolean isValid() {
		return isValid;
	}
}
