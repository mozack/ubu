package edu.unc.bioinf.ubu.fastq;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Creates fastq files in a format acceptable to Mapsplice.
 * Mapsplice currently does not accept fastq ids that contain spaces or do not end with /1 or /2.
 * Casava 1.8 is currently outputting fastq files in the above format.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqFormatter {

    private FastqInputFile input;
    private FastqOutputFile output;
    private String idSuffix;
    private boolean shouldConvertPhred33To64;
    private boolean shouldStripAfterSpace;
    
    public FastqFormatter(FastqFormatterOptions options) throws FileNotFoundException, IOException {    	
    	if (options.hasSuffix()) {
    		this.idSuffix     = options.getSuffix();
    	}
    	
    	this.shouldConvertPhred33To64 = options.shouldConvertPhred33To64();
    	this.shouldStripAfterSpace = options.shouldStripAfterWhitespace();
    	input.init(options.getInputFile());
    	output.init(options.getOutputFile());
    }

    FastqFormatter(FastqInputFile input, FastqOutputFile output, String idSuffix, boolean shouldStripAfterSpace, boolean shouldConvertPhred33To64) {
        this.input = input;
        this.output = output;
        this.idSuffix = idSuffix;
        this.shouldConvertPhred33To64 = shouldConvertPhred33To64;
        this.shouldStripAfterSpace = shouldStripAfterSpace;
    }
    
    public void process() throws IOException {
        FastqRecord rec = input.getNextRecord();
        
        int count = 0;
        
        while (rec != null) {
        	if (shouldStripAfterSpace) {
        		rec.stripNonReadInfoInId();
        	}
        	
        	if (idSuffix != null) {
        		rec.appendToId(idSuffix);
        	}
        	
        	if (shouldConvertPhred33To64) {
        		rec.phred33To64();
        	}
        	
            output.write(rec);
            rec = input.getNextRecord();
            
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " records.");
            }
        }
        
        input.close();
        output.close();
        
        System.out.println("Done.");
    }
    
    public static void run(String[] args) throws IOException {
    	FastqFormatterOptions options = new FastqFormatterOptions();
    	options.parseOptions(args);
        
    	if (options.isValid()) {
	        FastqFormatter prep = new FastqFormatter(options);
	        
	        prep.process();
    	}
    }    
}
