package edu.unc.bioinf.ubu.fastq;

import joptsimple.OptionParser;
import edu.unc.bioinf.ubu.util.Options;

/**
 * Options parser for {@code FastqMapsplicePrep}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqMapsplicePrepOptions extends Options {
	
	private static final String INPUT = "in";
	private static final String OUTPUT = "out";
	private static final String SUFFIX = "suffix";

	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT, "Input FASTQ file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT, "Output FASTQ file").withRequiredArg().ofType(String.class);
            parser.accepts(SUFFIX, "Read suffix (i.e. /1 or /2)").withRequiredArg().ofType(String.class);
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
        
        if (!getOptions().hasArgument(OUTPUT)) {
            isValid = false;
            System.err.println("Missing required output FASTQ file");
        }
        
        if (!getOptions().hasArgument(SUFFIX)) {
            isValid = false;
            System.err.println("Missing required read suffix");
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

	public String getSuffix() {
		return (String) getOptions().valueOf(SUFFIX);
	}
	
    public boolean isValid() {
        return isValid;
    }
}
