package edu.unc.bioinf.ubu.fastq;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code Sam2Fastq}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2FastqOptions extends Options {
    private static final String INPUT_SAM = "in";
    private static final String FASTQ1 = "fastq1";
    private static final String FASTQ2 = "fastq2";
    private static final String END1_SUFFIX = "end1";
    private static final String END2_SUFFIX = "end2";
    private static final String MAPSPLICE_FUSIONS = "mapsplice";
    
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(FASTQ1, "Required output FASTQ file").withRequiredArg().ofType(String.class);
            parser.accepts(FASTQ2, "Second FASTQ file for paired end").withRequiredArg().ofType(String.class);
            parser.accepts(END1_SUFFIX, "Id suffix used to identify the first read in a pair.  i.e. /1 (omit this option to use bit flag)").withRequiredArg().ofType(String.class);
            parser.accepts(END2_SUFFIX, "Id suffix used to identify the second read in a pair.  i.e. /2 (omit this option to use bit flag)").withRequiredArg().ofType(String.class);
            parser.accepts(MAPSPLICE_FUSIONS, "Enables special handling of Mapsplice fusions");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT_SAM)) {
            isValid = false;
            System.err.println("Missing required input SAM/BAM file");
        }
        
        if (!getOptions().hasArgument(FASTQ1)) {
            isValid = false;
            System.err.println("Missing required output SAM/BAM file");
        }
        
        //
        //  Identifying reads by name only makes sense for paired end.
        if ((getOptions().has(END1_SUFFIX)) &&
        	(!getOptions().has(FASTQ2))) {
        	isValid = false;
        	System.err.println(END1_SUFFIX + " only applicable for paired end.");
        }
        
        //
        //  Identifying reads by name only makes sense for paired end.
        if ((getOptions().has(END2_SUFFIX)) &&
        	(!getOptions().has(FASTQ2))) {
        	isValid = false;
        	System.err.println(END2_SUFFIX + " only applicable for paired end.");
        }
        
        //  Either specify both suffixes or neither.
        if ((getOptions().has(END1_SUFFIX) && !getOptions().has(END2_SUFFIX)) ||
        	(!getOptions().has(END1_SUFFIX) && getOptions().has(END2_SUFFIX))) {
        	isValid = false;
        	System.err.println("Please either specify both " + END1_SUFFIX + " and " + END2_SUFFIX + " or neither.");
        }
        
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT_SAM);
	}
	
	public String getFastq1() {
		return (String) getOptions().valueOf(FASTQ1);
	}
	
	public String getFastq2() {
		return (String) getOptions().valueOf(FASTQ2);
	}
	
	public boolean isPairedEnd() {
		return getOptions().hasArgument(FASTQ2);
	}
	
	public boolean shouldIdEndByReadName() {
		return getOptions().hasArgument(END1_SUFFIX);
	}
	
	public boolean isMapspliceFusions() {
		return getOptions().has(MAPSPLICE_FUSIONS);
	}
	
	public String getEnd1Suffix() {
		return (String) getOptions().valueOf(END1_SUFFIX);
	}
	
	public String getEnd2Suffix() {
		return (String) getOptions().valueOf(END2_SUFFIX);
	}
		
    public boolean isValid() {
        return isValid;
    }
}