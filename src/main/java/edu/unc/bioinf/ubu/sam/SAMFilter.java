package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Filters paired reads from a SAM or BAM file.  Indels and/or clusters
 * greater than a specified insert length may be filtered.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SAMFilter {
    
    private boolean shouldStripIndels = false;
    private int     maxInsertLen = -1;
    private int     minMappingQuality = -1;

    public void filter(String input, String output) {
        File outputFile = new File(output);
        SamReadPairReader reader = new SamReadPairReader(input);
        
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getHeader(),
                false, outputFile);
        
        for (ReadPair pair : reader) {
            SAMRecord read1 = pair.getRead1();
            SAMRecord read2 = pair.getRead2();
            
            boolean isPairIncluded = !shouldStripIndels || !hasIndels(read1, read2);
            
            if (isPairIncluded) {
                isPairIncluded = !isMaxInsertLenExceeded(read1, read2);
            }
            
            if (isPairIncluded) {
            	isPairIncluded = !isBelowMinMappingQuality(read1, read2);
            }
            
            // Output only the pairs that have passed our tests
            if (isPairIncluded) {
                writer.addAlignment(read1);
                writer.addAlignment(read2);
            }
        }
        
        writer.close();
        reader.close();
    }
    
    private boolean isBelowMinMappingQuality(SAMRecord read1, SAMRecord read2) {
    	boolean isBelowMin = false;
    	
    	if (isMinMappingQualitySpecified()) {
    		isBelowMin = 
    			(read1.getMappingQuality() < minMappingQuality) ||
    			(read2.getMappingQuality() < minMappingQuality);
    	}
    	
    	return isBelowMin;
    }
    
    private boolean isMaxInsertLenExceeded(SAMRecord read1, SAMRecord read2) {
        boolean isMaxExceeded = false;
        
        if (isMaxInsertLenSpecified()) {
            isMaxExceeded = 
                (Math.abs(read1.getInferredInsertSize()) > maxInsertLen) &&
                (Math.abs(read2.getInferredInsertSize()) > maxInsertLen);
        }
        
        return isMaxExceeded;
    }
    
    private void setMaxInsertLen(int len) {
        maxInsertLen = len;
    }
    
    private boolean isMaxInsertLenSpecified() {
        return maxInsertLen > 0;
    }
    
    private void setMinMappingQuality(int min) {
    	minMappingQuality = min;
    }
    
    private boolean isMinMappingQualitySpecified() {
    	return minMappingQuality > 0;
    }
    
    private boolean hasIndels(SAMRecord read1, SAMRecord read2) {
        return ( (hasIndel(read1)) || (hasIndel(read2)) );
    }
    
    public void setShouldStripIndels(boolean shouldStripIndels) {
        this.shouldStripIndels = true;
    }
    
    private boolean hasIndel(SAMRecord read) {
        for (CigarElement element : read.getCigar().getCigarElements()) {
            if ((element.getOperator() == CigarOperator.D) ||
                (element.getOperator() == CigarOperator.I)) {
                return true;
            }
        }
        
        return false;
    }
        
    public static void run(String[] args) throws IOException {
    	SamFilterOptions options = new SamFilterOptions();
    	options.parseOptions(args);
    	
    	if (options.isValid()) {
    		long s = System.currentTimeMillis();
    		
    		SAMFilter filter = new SAMFilter();
    		
    		filter.setMaxInsertLen(options.getMaxInsertLen());
    		filter.setMinMappingQuality(options.getMinMappingQuality());
    		filter.setShouldStripIndels(options.shouldStripIndels());
    		filter.filter(options.getInputFile(), options.getOutputFile());
    		
            long e = System.currentTimeMillis();
            
            System.out.println("Elapsed: " + (e-s)/1000);
    	}
    }
}
