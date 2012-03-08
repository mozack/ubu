package edu.unc.bioinf.ubu.fastq;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Creates fastq files in a format acceptable to Mapsplice.
 * Mapsplice currently does not accept fastq ids that contain spaces or do not end with /1 or /2.
 * Casava 1.8 is currently outputting fastq files in the above format.
 * 
 * @author lmose
 */
public class FastqMapsplicePrep {

    private FastqInputFile input;
    private FastqOutputFile output;
    private String idSuffix;
    
    public FastqMapsplicePrep(String inputFile, String outputfile, String idSuffix) throws FileNotFoundException, IOException {
        this(new FastqInputFile(), new FastqOutputFile(), idSuffix); 
        input.init(inputFile);
        output.init(outputfile);
    }
    
    FastqMapsplicePrep(FastqInputFile input, FastqOutputFile output, String idSuffix) {
        this.input = input;
        this.output = output;
        this.idSuffix = idSuffix;
    }
    
    public void process() throws IOException {
        FastqRecord rec = input.getNextRecord();
        
        int count = 0;
        
        while (rec != null) {
            rec.stripNonReadInfoInId();
            rec.appendToId(idSuffix);
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
    	FastqMapsplicePrepOptions options = new FastqMapsplicePrepOptions();
    	options.parseOptions(args);
        
    	if (options.isValid()) {
	        FastqMapsplicePrep prep = new FastqMapsplicePrep(options.getInputFile(), options.getOutputFile(),
	        		options.getSuffix());
	        
	        prep.process();
    	}
    }
    
    /*
    public static void main(String[] args) throws IOException {
//        String input = "/home/lisle/mapsplice12/out.fastq";
//        String output = "/home/lisle/mapsplice12/out2.fastq";
//        String suffix = "/1";
        
//        if (args.length != 3) {
//            System.out.println("Usage: FastqMapsplicePrep <input_file> <output_file> <suffix>");
//            System.exit(-1);
//        }
//        
//        String input = args[0];
//        String output = args[1];
//        String suffix = args[2];
    	
    }
    */
}
