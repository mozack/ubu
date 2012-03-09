package edu.unc.bioinf.ubu.sam;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code SamSummarizer}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamSummarizerOptions extends Options {
	
	private static final String INPUT = "in";
	private static final String OUTPUT = "out";
	private static final String HEADER = "header";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT, "Input SAM/BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT, "Output summary file").withRequiredArg().ofType(String.class);
            parser.accepts(HEADER, "Output header");
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
            System.err.println("Missing required output summary file");
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
	
	public boolean shouldOutputHeader() {
		return getOptions().has(HEADER);
	}
	
    public boolean isValid() {
        return isValid;
    }
}
