package edu.unc.bioinf.ubu.sam;

import java.io.IOException;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * Handles parsing of command line arguments for the GenomeToTranscriptomeConverter
 * 
 * @author lmose
 */
public class GenomeToTranscriptomeOptions {
    
    // Option flags
    private static final String BED_FILE    = "bed";
    private static final String READ_OFFSET = "offset";
    private static final String INPUT_ALIGNMENT_FILE = "in";
    private static final String OUTPUT_ALIGNMENT_FILE = "out";
    private static final String ORDERING_FASTA = "order";
    private static final String REVERSE_STRAND_COORDS = "reverse";
    private static final String OUTPUT_XG_TAGS = "xgtags";
    private static final String SINGLE_END = "single";    
    private static final String HELP = "help";
    
    private static final int DEFAULT_READ_OFFSET = 25;
    
    private OptionSet options;
    private OptionParser parser;
    private boolean isValid = false;

    public void parseOptions(String[] args) {
        parser = new OptionParser();
        parser.accepts(BED_FILE, "Bed file definition of isoforms").withRequiredArg().ofType(String.class);
        parser.accepts(INPUT_ALIGNMENT_FILE, "Input alignment file").withRequiredArg().ofType(String.class);
        parser.accepts(OUTPUT_ALIGNMENT_FILE, "Output alignment file").withRequiredArg().ofType(String.class);
        parser.accepts(READ_OFFSET, "Optional Offset size of read index (default 25)").withRequiredArg().ofType(Integer.class);
        parser.accepts(ORDERING_FASTA, "Optional FASTA file used to determine order of isoforms in BAM header (important for RSEM)").withRequiredArg().ofType(String.class);
        parser.accepts(REVERSE_STRAND_COORDS, "Optional flag indicating that reverse strand coordinates should be reported");
        parser.accepts(OUTPUT_XG_TAGS, "Optional flag indicating that genomic coordinates should be output in a XG tag");
        parser.accepts(SINGLE_END, "Optional flag indicating that reads need not be paired in the same transcript to be output (default is off)");
        parser.accepts(HELP, "Print this help message");

        try {
            options = parser.parse(args);
        
            if (options.has(HELP)) {
                printHelp();
            } else {
                validate();
            }
        } catch (joptsimple.OptionException e) {
            System.err.println(e.getMessage());
            printHelp();
        }
    }
    
    private void printHelp() {
        try {
            parser.printHelpOn(System.err);
        }
        catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("IOException encountered when attempting to output help.");
        }
    }
    
    private void validate() {
        isValid = true;
        
        if (!options.hasArgument(BED_FILE)) {
            isValid = false;
            System.err.println("Missing required bed file");
        }
        
        if (!options.hasArgument(INPUT_ALIGNMENT_FILE)) {
            isValid = false;
            System.err.println("Missing required input alignment file");
        }
        
        if (!options.hasArgument(OUTPUT_ALIGNMENT_FILE)) {
            isValid = false;
            System.err.println("Missing required output alignment file");
        }
        
        if (!isValid) {
            printHelp();
        }
    }
    
    public boolean isValid() {
        return isValid;
    }
    
    public String getBedFile() {
        return (String) options.valueOf(BED_FILE);
    }
    
    public int getReadOffset() {
        int readOffset = DEFAULT_READ_OFFSET;
        
        if (options.hasArgument(READ_OFFSET)) {
            readOffset = (Integer) options.valueOf(READ_OFFSET);
        }
        
        return readOffset;
    }
    
    public String getInputAlignmentFile() {
        return (String) options.valueOf(INPUT_ALIGNMENT_FILE);
    }
    
    public String getOutputAlignmentFile() {
        return (String) options.valueOf(OUTPUT_ALIGNMENT_FILE);
    }
    
    public String getOrderingFastaFile() {
        return (String) options.valueOf(ORDERING_FASTA);
    }
    
    public boolean hasOrderingFastaFile() {
        return options.hasArgument(ORDERING_FASTA);
    }
    
    public boolean isPositiveStrandReportingOnly() {
        return !options.has(REVERSE_STRAND_COORDS);
    }
    
    public boolean shouldOutputXgTags() {
        return options.has(OUTPUT_XG_TAGS);
    }
    
    public boolean isSingleEnd() {
    	return options.has(SINGLE_END);
    }
}
