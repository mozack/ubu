package edu.unc.bioinf.ubu.fastq;

import java.io.IOException;

/**
 * Splits a single fastq file containing concatenated paired end reads
 * into 2 fastq files.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqSplitter {
	
	public void split(String input, String output1, String output2) throws IOException {
	    FastqInputFile in = new FastqInputFile();
	    FastqOutputFile out1 = new FastqOutputFile();
	    FastqOutputFile out2 = new FastqOutputFile();
	    
	    FastqRecord rec = in.getNextRecord();
	    
	    while (rec != null) {
	    	int len = rec.getSequence().length();
	    	
	    	String id1 = rec.getId() + "/1";
	    	String bases1 = rec.getSequence().substring(0, len/2);
	    	String qualities1 = rec.getQuality().substring(0, len/2);
	    	FastqRecord outRec1 = new FastqRecord(id1, bases1, qualities1);
	    	out1.write(outRec1);
	    	
	    	String id2 = rec.getId() + "/2";
	    	String bases2 = rec.getSequence().substring(len/2, len);
	    	String qualities2 = rec.getQuality().substring(len/2, len);
	    	FastqRecord outRec2 = new FastqRecord(id2, bases2, qualities2);
	    	out2.write(outRec2);
	    }
	    
	    in.close();
	    out1.close();
	    out2.close();
	}
	
	public static void run(String[] args) throws IOException {
		FastqSplitterOptions options = new FastqSplitterOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			FastqSplitter splitter = new FastqSplitter();
		    
			splitter.split(options.getInputFile(), options.getOutputFile1(), options.getOutputFile2());
		}
	}
}
