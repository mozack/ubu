package edu.unc.bioinf.ubu.sam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Counts splice junctions in a bam or sam file
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SpliceJunctionCounter {
    
    private SpliceJunctionMap spliceJunctionMap;
    private Map<SpliceJunction, Integer> spliceJunctionCounts = new HashMap<SpliceJunction, Integer>(); 
    
    public SpliceJunctionCounter(SpliceJunctionMap spliceJunctionMap) {
        this.spliceJunctionMap = spliceJunctionMap;
    }
    
    public void count(String inputFile, String outputFile) throws IOException {
        
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
        
        File file = new File(inputFile);
        
        SAMFileReader inputSam = new SAMFileReader(file);
        inputSam.setValidationStringency(ValidationStringency.SILENT);
        
        int count = 0;

        for (SAMRecord read : inputSam) {
            updateJunctionCounts(read);
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
        }
                
        outputCounts(writer);
        writer.close();
    }
    
    private void outputCounts(BufferedWriter writer) throws IOException {
        System.out.println("Writing counts.");
        
        for (SpliceJunction junction : spliceJunctionMap.getAllSpliceJunctions()) {
            String junctionKey = spliceJunctionMap.getJunctionKey(junction);
            int count = 0;
            
            if (spliceJunctionCounts.containsKey(junction)) {
                count = spliceJunctionCounts.get(junction);
            }
            
            String line = junctionKey + "\t" + count + "\n";

            writer.write(line);
        }
    }
    
    private void updateCount(SpliceJunction junction) {
        Integer count = spliceJunctionCounts.get(junction);
        if (count == null) {
            spliceJunctionCounts.put(junction, 1);
        } else {
            spliceJunctionCounts.put(junction, count+1);
        }
    }
    
    private void updateJunctionCounts(SAMRecord read) {
        
        if (read.getCigarString().contains("N")) {
            for (ReadBlock block : ReadBlock.getReadBlocks(read)) {
                if (block.getType() == CigarOperator.N) {
                    updateCount(new SpliceJunction(read.getReferenceName(), block.getReferenceStart(), block.getReferenceStop()));
                }
            }
        }
    }
    
    public static void run(String[] args) throws IOException {
    	SpliceJunctionCounterOptions options = new SpliceJunctionCounterOptions();
    	options.parseOptions(args);
    	
    	if (options.isValid()) {
    		
            long start = System.currentTimeMillis();
            
    		SpliceJunctionMap map = new SpliceJunctionMap(options.getJunctionFile());
    		SpliceJunctionCounter counter = new SpliceJunctionCounter(map);
    		counter.count(options.getInputFile(), options.getOutputFile());
    		
            long stop = System.currentTimeMillis();
            
            System.out.println("Done.  Elapsed secs: " + (stop-start)/1000);
    	}
    }
    
    public static void main(String[] args) throws Exception {
        SpliceJunctionMap map = new SpliceJunctionMap("/home/lisle/gaf/splice_junctions.txt");
        SpliceJunctionCounter counter = new SpliceJunctionCounter(map);
        
        counter.count("/home/lisle/data/junction/small_sorted_by_read.bam", "/home/lisle/data/junction/counts2.txt"); 
//        counter.count("/home/lisle/data/junction/1.sam", "/home/lisle/data/junction/counts.txt");
        
        
//        SpliceJunctionMap map = new SpliceJunctionMap(args[0]);
//        SpliceJunctionCounter counter = new SpliceJunctionCounter(map);
//        
//        counter.count(args[1], args[2]); 
    }
}
