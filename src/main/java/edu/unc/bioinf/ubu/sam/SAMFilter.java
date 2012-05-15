package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Filters reads from a SAM or BAM file.  For paired end, read pairs are filtered.
 * Candidates for filtering include:<br/>
 * Indels<br/>
 * Clusters greater than a specified insert length<br/>
 * Reads with low mapping quality<br/>
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SAMFilter {
    
	private boolean isPairedEnd = true;
    private boolean shouldStripIndels = false;
    private int     maxInsertLen = -1;
    private int     minMappingQuality = -1;

    public void filter(String input, String output) {
    	if (isPairedEnd) {
    		filterPairedEnd(input, output);
    	} else {
    		filterSingleEnd(input, output);
    	}
    }
    
    private void filterPairedEnd(String input, String output) {
        File outputFile = new File(output);
        SamReadPairReader reader = new SamReadPairReader(input);
        
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getHeader(),
                true, outputFile);
        
        for (ReadPair pair : reader) {
            // Output only the pairs that have passed our tests
            if (isReadPairIncluded(pair)) {
                writer.addAlignment(pair.getRead1());
                writer.addAlignment(pair.getRead2());
            }
        }
        
        writer.close();
        reader.close();
    }
    
    private void filterSingleEnd(String input, String output) {
        File outputFile = new File(output);
        SAMFileReader reader = new SAMFileReader(new File(input));
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(),
                true, outputFile);
        
        int cnt = 0;
        for (SAMRecord read : reader) {
            // Output only the pairs that have passed our tests
            if (isReadIncluded(read)) {
                writer.addAlignment(read);
            }
            cnt++;
            if ((cnt % 1000000) == 0) {
                System.out.println("Processed reads: " + cnt);
            }
        }
        
        writer.close();
        reader.close();
    }
    
    boolean isReadPairIncluded(ReadPair pair) {
    	if (!isPairedEnd) {
    		throw new UnsupportedOperationException("Invalid call to isReadPairIncluded for single end filtering.");
    	}
    	
    	return isReadIncluded(pair.getRead1()) && isReadIncluded(pair.getRead2());
    }
    
    boolean isReadIncluded(SAMRecord read) {
    	return
    		((!hasIndel(read)) &&
    		 (!isMaxInsertLenExceeded(read)) &&
    		 (!isBelowMinMappingQuality(read)));
    }
    
    private boolean isBelowMinMappingQuality(SAMRecord read) {
    	boolean isBelowMin = false;
    	
    	if (isMinMappingQualitySpecified()) {
    		isBelowMin = (read.getMappingQuality() < minMappingQuality);
    	}
    	
    	return isBelowMin;
    }
    
    private boolean isMaxInsertLenExceeded(SAMRecord read) {
        boolean isMaxExceeded = false;
        
        if (isMaxInsertLenSpecified()) {
            isMaxExceeded = (Math.abs(read.getInferredInsertSize()) > maxInsertLen);
        }
        
        return isMaxExceeded;
    }
    
    public void setMaxInsertLen(int len) {
        maxInsertLen = len;
    }
    
    private boolean isMaxInsertLenSpecified() {
        return maxInsertLen > 0;
    }
    
    public void setMinMappingQuality(int min) {
    	minMappingQuality = min;
    }
    
    private boolean isMinMappingQualitySpecified() {
    	return minMappingQuality > 0;
    }
    
    public void setShouldStripIndels(boolean shouldStripIndels) {
        this.shouldStripIndels = true;
    }
    
    public void setPairedEnd(boolean isPairedEnd) {
		this.isPairedEnd = isPairedEnd;
	}

	private boolean hasIndel(SAMRecord read) {
    	if (shouldStripIndels) {
	        for (CigarElement element : read.getCigar().getCigarElements()) {
	            if ((element.getOperator() == CigarOperator.D) ||
	                (element.getOperator() == CigarOperator.I)) {
	                return true;
	            }
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
    		
    		filter.setPairedEnd(options.isPairedEnd());
    		filter.setMaxInsertLen(options.getMaxInsertLen());
    		filter.setMinMappingQuality(options.getMinMappingQuality());
    		filter.setShouldStripIndels(options.shouldStripIndels());
    		filter.filter(options.getInputFile(), options.getOutputFile());
    		
            long e = System.currentTimeMillis();
            
            System.out.println("Elapsed: " + (e-s)/1000);
    	}
    }
}
