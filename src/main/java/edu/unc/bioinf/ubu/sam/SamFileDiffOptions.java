package edu.unc.bioinf.ubu.sam;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code SamFileDiff}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamFileDiffOptions extends Options {
	
	private static final String INPUT1  = "in1";
	private static final String INPUT2  = "in2";
	private static final String OUTPUT1 = "out1";
	private static final String OUTPUT2 = "out2";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT1, "Input SAM/BAM file 1").withRequiredArg().ofType(String.class);
            parser.accepts(INPUT2, "Input SAM/BAM file 2").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT1, "Output SAM/BAM file containing reads unique to input file 1").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT2, "Output SAM/BAM file containing reads unique to input file 2").withRequiredArg().ofType(String.class);
            parser.accepts(HELP, "Print this help message");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT1)) {
            isValid = false;
            System.err.println("Missing required input file 1");
        }
        
        if (!getOptions().hasArgument(INPUT2)) {
            isValid = false;
            System.err.println("Missing required input file 2");
        }
        
        if (!getOptions().hasArgument(OUTPUT1)) {
            isValid = false;
            System.err.println("Missing required output file 1");
        }
        
        if (!getOptions().hasArgument(OUTPUT2)) {
            isValid = false;
            System.err.println("Missing required output file 2");
        }
        
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInput1File() {
		return (String) getOptions().valueOf(INPUT1);
	}
	
	public String getInput2File() {
		return (String) getOptions().valueOf(INPUT2);
	}
	
	public String getOutput1File() {
		return (String) getOptions().valueOf(OUTPUT1);
	}
	
	public String getOutput2File() {
		return (String) getOptions().valueOf(OUTPUT2);
	}
	
    public boolean isValid() {
        return isValid;
    }
}
