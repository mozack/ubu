package edu.unc.bioinf.ubu.fastq;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code FastqMapsplicePrep}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqSplitterOptions extends Options {
	
	private static final String INPUT = "in";
	private static final String OUTPUT1 = "out1";
	private static final String OUTPUT2 = "out2";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT, "Input FASTQ file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT1, "Output FASTQ file 1").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT2, "Output FASTQ file 2").withRequiredArg().ofType(String.class);
            parser.accepts(HELP, "Print this help message");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
        isValid = true;
        
        if (!getOptions().hasArgument(INPUT)) {
            isValid = false;
            System.err.println("Missing required input FASTQ file");
        }
        
        if (!getOptions().hasArgument(OUTPUT1)) {
            isValid = false;
            System.err.println("Missing required output FASTQ file 1");
        }
        
        if (!getOptions().hasArgument(OUTPUT2)) {
            isValid = false;
            System.err.println("Missing required output FASTQ file 1");
        }
        
        if (!isValid) {
            printHelp();
        }
	}
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT);
	}
	
	public String getOutputFile1() {
		return (String) getOptions().valueOf(OUTPUT1);
	}
	
	public String getOutputFile2() {
		return (String) getOptions().valueOf(OUTPUT2);
	}
	
    public boolean isValid() {
        return isValid;
    }
}
