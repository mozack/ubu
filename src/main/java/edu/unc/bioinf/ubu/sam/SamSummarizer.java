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
 * Provides summary statistics per reference for a SAM/BAM file.  Unmapped
 * reads are not included in summary stats.
 * 
 * Output columns:
 * 
 * reference name
 * number of aligned bases
 * edit distance count
 * base error rate (<edit distance count>/<number of aligned bases>)
 * number of reads
 * number of reads with zero mapping quality
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamSummarizer {
	
	private static final String EDIT_DISTANCE_TAG = "NM";
	
	private Map<String, ReferenceCounts> refCountMap = new HashMap<String, ReferenceCounts>();
	
    public void summarize(String inputFile, String outputFile, boolean shouldOutputHeader) throws IOException {
        
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

        if (shouldOutputHeader) {
        	outputHeader(writer);
        }
        
        outputAllCounts(writer);
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
	    	
	    	if (read.getMappingQuality() == 0) {
	    		counts.incrementMappingQualityZeroCount(1);
	    	}
	    	
	    	counts.incrementReadCount(1);
    	}
    }
    
    private void outputAllCounts(BufferedWriter writer) throws IOException {
    	
    	ReferenceCounts totals = new ReferenceCounts();
    	
    	for (String ref : getSortedReferences()) {
    		ReferenceCounts counts = getReferenceCounts(ref);
    		outputCounts(counts, ref, writer);
    		
    		totals.incrementAlignedBases(counts.getAlignedBases());
    		totals.incrementEditDistanceCount(counts.getEditDistanceCount());
    		totals.incrementMappingQualityZeroCount(counts.getMappingQualityZeroCount());
    		totals.incrementReadCount(counts.getReadCount());
    	}
    	
    	outputCounts(totals, "Total", writer);
    }
    
    private void outputHeader(BufferedWriter writer) throws IOException {
    	StringBuffer buf = new StringBuffer();
    	
		buf.append("ref");
		buf.append('\t');
		buf.append("Aligned_Bases");
		buf.append('\t');
		buf.append("NM");
		buf.append('\t');
		buf.append("Error_Rate");
		buf.append('\t');
		buf.append("Reads");
		buf.append('\t');
		buf.append("0_Mapping_Quality");
		buf.append('\t');
		buf.append("0_Mapping_Quality_Rate");
		buf.append('\n');

		writer.write(buf.toString());
    }
    
    private void outputCounts(ReferenceCounts counts, String ref, BufferedWriter writer) throws IOException {
		StringBuffer buf = new StringBuffer();
		buf.append(ref);
		buf.append('\t');
		buf.append(counts.getAlignedBases());
		buf.append('\t');
		buf.append(counts.getEditDistanceCount());
		buf.append('\t');
		buf.append(counts.getErrorRate());
		buf.append('\t');
		buf.append(counts.getReadCount());
		buf.append('\t');
		buf.append(counts.getMappingQualityZeroCount());
		buf.append('\t');
		buf.append(counts.getMappingQualityZeroRate());
		buf.append('\n');
		
		writer.write(buf.toString());
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
    	private long mappingQualityZeroCount = 0;
    	private long readCount = 0;
    	
    	void incrementAlignedBases(long inc) {
    		alignedBases += inc;
    	}
    	
    	void incrementEditDistanceCount(long inc) {
    		editDistanceCount += inc;
    	}
    	
    	void incrementMappingQualityZeroCount(long inc) {
    		mappingQualityZeroCount += inc;
    	}
    	
    	void incrementReadCount(long inc) {
    		readCount += inc;
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
    	
    	long getMappingQualityZeroCount() {
    		return mappingQualityZeroCount;
    	}
    	
    	long getReadCount() {
    		return readCount;
    	}
    	
    	double getMappingQualityZeroRate() {
    		return (double) mappingQualityZeroCount / (double) readCount;
    	}
    }
    
    public static void run(String[] args) throws IOException {
    	SamSummarizerOptions options = new SamSummarizerOptions();
    	options.parseOptions(args);
    	
    	if (options.isValid()) {
    		new SamSummarizer().summarize(options.getInputFile(),
    				options.getOutputFile(), options.shouldOutputHeader());
    	}
    }
    
    public static void main(String[] args) throws IOException {
    	run(args);    	
    }
}
