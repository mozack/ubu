package edu.unc.bioinf.ubu.sam;

import edu.unc.bioinf.ubu.util.Options;
import joptsimple.OptionParser;

/**
 * Options for {@code SamConverter}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamConverterOptions extends Options {
	private static final String INPUT  = "in";
	private static final String OUTPUT = "out";
	private static final String PHRED64_TO_PHRED33 = "phred64to33";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT, "Input SAM/BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT, "Output SAM/BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(PHRED64_TO_PHRED33, "If specified, convert quality score from phred64 to phred33");
            parser.accepts(HELP, "Print this help message");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT)) {
            isValid = false;
            System.err.println("Missing required input SAM/BAM file");
        }
        
        if (!getOptions().hasArgument(OUTPUT)) {
            isValid = false;
            System.err.println("Missing required output SAM/BAM file");
        }
        
        if (!getOptions().has(PHRED64_TO_PHRED33)) {
            isValid = false;
            System.err.println("No conversion options specified.");
        }
                
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT);
	}
	
	public String getOutputFile() {
		return (String) getOptions().valueOf(OUTPUT);
	}

	public boolean shouldConvertPhred64toPhred33() {
		return getOptions().has(PHRED64_TO_PHRED33);
	}
		
    public boolean isValid() {
        return isValid;
    }
}
