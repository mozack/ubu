package edu.unc.bioinf.ubu.sam;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code SpliceJunctionCounter}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SpliceJunctionCounterOptions extends Options {
	
	private static final String INPUT  = "in";
	private static final String OUTPUT = "out";
	private static final String JUNCTION_FILE = "junctions";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT, "Input SAM/BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT, "Output file containing junction counts").withRequiredArg().ofType(String.class);
            parser.accepts(JUNCTION_FILE, "Input list of junctions defined in a format similar to: chr1:12227:+,chr1:12595:+").withRequiredArg().ofType(String.class);
            parser.accepts(HELP, "Print this help message");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT)) {
            isValid = false;
            System.err.println("Missing required input file 1");
        }
        
        if (!getOptions().hasArgument(OUTPUT)) {
            isValid = false;
            System.err.println("Missing required input file 2");
        }
        
        if (!getOptions().hasArgument(JUNCTION_FILE)) {
            isValid = false;
            System.err.println("Missing required input junction file");
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
	
	public String getJunctionFile() {
		return (String) getOptions().valueOf(JUNCTION_FILE);
	}
		
    public boolean isValid() {
        return isValid;
    }
}
