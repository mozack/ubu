package edu.unc.bioinf.ubu.sam;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

public class SamFilterOptions extends Options {
    private static final String INPUT_FILE = "in";
    private static final String OUTPUT_FILE = "out";
    private static final String STRIP_INDELS = "strip-indels";
    private static final String MAX_INSERT_LEN = "max-insert";
    private static final String MAPPING_QUALITY = "mapq";
    private static final String SINGLE_END = "single";
    private static final String INCLUDE_INDELS_ONLY = "indels-only";
    
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_FILE, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_FILE, "Required output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(SINGLE_END, "If specified, process bam as single end, discarding reads independently (default paired end)");
            parser.accepts(STRIP_INDELS, "If specified, discard read pairs containing indels from output (default off)");
            parser.accepts(MAX_INSERT_LEN, "If specified, discard clusters greater than specified insert length").withRequiredArg().ofType(Integer.class);
            parser.accepts(MAPPING_QUALITY, "If specified, discard clusters with mapping quality less than the specified value").withRequiredArg().ofType(Integer.class);
            parser.accepts(INCLUDE_INDELS_ONLY, "If specified, discard reads not containing indels (default off)");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT_FILE)) {
            isValid = false;
            System.err.println("Missing required input SAM/BAM file");
        }
        
        if (!getOptions().hasArgument(OUTPUT_FILE)) {
            isValid = false;
            System.err.println("Missing required output SAM/BAM file");
        }
        
        // Validate that there is filtering to be done.
        // If no options are specified, only paired reads will be output.
        if ((!getOptions().has(STRIP_INDELS)) && 
        	(!getOptions().hasArgument(MAX_INSERT_LEN)) &&
        	(!getOptions().hasArgument(MAPPING_QUALITY)) &&
        	(!getOptions().has(INCLUDE_INDELS_ONLY)) &&
        	(getOptions().has(SINGLE_END))) {
        	isValid = false;
        	System.err.println("At least one filtering option must be specified");
        }
        
        if ((getOptions().has(STRIP_INDELS)) && (getOptions().has(INCLUDE_INDELS_ONLY))) {
        	isValid = false;
        	System.err.println("Cannot specify both " + STRIP_INDELS + " and " + INCLUDE_INDELS_ONLY);
        }
        
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT_FILE);
	}
	
	public String getOutputFile() {
		return (String) getOptions().valueOf(OUTPUT_FILE);
	}
	
	public boolean shouldStripIndels() {
		return getOptions().has(STRIP_INDELS);
	}
	
	public boolean shouldIncludeIndelsOnly() {
		return getOptions().has(INCLUDE_INDELS_ONLY);
	}
	
	public boolean isPairedEnd() {
		return !getOptions().has(SINGLE_END);
	}
	
	public int getMaxInsertLen() {
		int maxLen = -1;
		
		if (getOptions().hasArgument(MAX_INSERT_LEN)) {
			maxLen = (Integer) getOptions().valueOf(MAX_INSERT_LEN);
		}
		
		return maxLen;
	}
	
	public int getMinMappingQuality() {
		int minMapQ = -1;
		
		if (getOptions().hasArgument(MAPPING_QUALITY)) {
			minMapQ = (Integer) getOptions().valueOf(MAPPING_QUALITY);
		}

		return minMapQ;
	}
	
    public boolean isValid() {
        return isValid;
    }
}
