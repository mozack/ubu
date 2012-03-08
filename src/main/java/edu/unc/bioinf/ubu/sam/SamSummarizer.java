package edu.unc.bioinf.ubu.sam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Summarizes aligned base and tag information per reference in a SAM or BAM file.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamSummarizer {
	
	private static final String EDIT_DISTANCE_TAG = "NM";
	
	private Map<String, ReferenceCounts> refCountMap = new HashMap<String, ReferenceCounts>();
	
    public void summarize(String inputFile, String outputFile) throws IOException {
        
        long start = System.currentTimeMillis();
        
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
        
        File file = new File(inputFile);
        
        SAMFileReader inputSam = new SAMFileReader(file);
        inputSam.setValidationStringency(ValidationStringency.SILENT);
        
        int count = 0;

        for (SAMRecord read : inputSam) {
        	
        	countRead(read);
        	
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
        }
                
        outputCounts(writer);
        writer.close();
        
        long stop = System.currentTimeMillis();
        
        System.out.println("free mem: " + Runtime.getRuntime().freeMemory());
        System.out.println("total mem: " + Runtime.getRuntime().totalMemory());
        
        System.out.println("Done.  Elapsed secs: " + (stop-start)/1000);
    }
    
    private void countRead(SAMRecord read) {
    	int alignedBases = 0;
    	if (!read.getReadUnmappedFlag()) {
    		for (ReadBlock block : ReadBlock.getReadBlocks(read)) {
    			if (block.getType() == CigarOperator.MATCH_OR_MISMATCH) {
    				alignedBases += block.getLength();
    			}
    		}
    	
	    	Integer editDistance = read.getIntegerAttribute(EDIT_DISTANCE_TAG);
	    	
	    	ReferenceCounts counts = getReferenceCounts(read.getReferenceName());
	    	counts.incrementAlignedBases(alignedBases);
	    	if (editDistance != null) {
	    		counts.incrementEditDistanceCount(editDistance);
	    	}
    	}
    }
    
    private void outputCounts(BufferedWriter writer) throws IOException {
    	for (String ref : getSortedReferences()) {
    		ReferenceCounts counts = getReferenceCounts(ref);
    		
    		StringBuffer buf = new StringBuffer();
    		buf.append(ref);
    		buf.append('\t');
    		buf.append(counts.getAlignedBases());
    		buf.append('\t');
    		buf.append(counts.getEditDistanceCount());
    		buf.append('\t');
    		buf.append(counts.getErrorRate());
    		buf.append('\n');
    		
    		writer.write(buf.toString());
    	}
    }
    
    
    private List<String> getSortedReferences() {
    	List<String> refs = new ArrayList<String>();
    	
    	refs.addAll(refCountMap.keySet());
    	Collections.sort(refs);
    	
    	return Collections.unmodifiableList(refs);
    }
    
    private ReferenceCounts getReferenceCounts(String reference) {
    	ReferenceCounts counts = refCountMap.get(reference);
    	
    	if (counts == null) {
    		counts = new ReferenceCounts();
    		refCountMap.put(reference, counts);
    	}
    	
    	return counts;
    }

    static class ReferenceCounts {
    	private long alignedBases = 0;
    	private long editDistanceCount = 0;
    	
    	void incrementAlignedBases(int inc) {
    		alignedBases += inc;
    	}
    	
    	void incrementEditDistanceCount(int inc) {
    		editDistanceCount += inc;
    	}
    	
    	long getAlignedBases() {
    		return alignedBases;
    	}
    	
    	long getEditDistanceCount() {
    		return editDistanceCount;
    	}
    	
    	double getErrorRate() {
    		return (double) editDistanceCount / (double) alignedBases;
    	}
    }
    
    public static void main(String[] args) throws IOException {
    	String input = args[0];
    	String output = args[1];
    	
//    	String input = "/home/lisle/data/summarizer/bwa1.sam";
//    	String output = "/home/lisle/data/summarizer/bwa1.tsv";
    	
    	new SamSummarizer().summarize(input, output);
    }
}
