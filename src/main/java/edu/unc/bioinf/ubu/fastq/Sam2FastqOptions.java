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
    private static final String ID_END_BY_READ_NAME = "use-name";
    
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(FASTQ1, "Required output FASTQ file").withRequiredArg().ofType(String.class);
            parser.accepts(FASTQ2, "Second FASTQ file for paired end").withRequiredArg().ofType(String.class);
            parser.accepts(ID_END_BY_READ_NAME, "If specified, use read name to identify end in pair.  i.e. /1 = first in pair (default off)");
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
        if ((getOptions().has(ID_END_BY_READ_NAME)) &&
        	(!getOptions().has(FASTQ2))) {
        	isValid = false;
        	System.err.println(ID_END_BY_READ_NAME + " only applicable for paired end.");
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
		return getOptions().has(ID_END_BY_READ_NAME);
	}
		
    public boolean isValid() {
        return isValid;
    }
}
