package edu.unc.bioinf.ubu.fastq;

import java.io.File;
import java.io.IOException;

import edu.unc.bioinf.ubu.sam.ReverseComplementor;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Converts SAM/BAM file to FASTQ
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2Fastq {
	
	private FastqOutputFile output1;
	private FastqOutputFile output2;
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	private boolean shouldIdentifyEndByReadId = false;

	/**
	 * Convert the input paired end SAM/BAM file into 2 fastq files.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFastq1, String outputFastq2) throws IOException {
		String last1Read = "";
		String last2Read = "";
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        output1 = new FastqOutputFile();
        output2 = new FastqOutputFile();
        
        int output1Count = 0;
        int output2Count = 0;
        
        for (SAMRecord read : reader) {
        	if (isFirstInPair(read)) {
        		if (!read.getReadName().equals(last1Read)) {
        			output1.write(samReadToFastqRecord(read));
        			last1Read = read.getReadName();
        			output1Count += 1;
        		}
        	} else {
        		if (!read.getReadName().equals(last2Read)) {
        			output2.write(samReadToFastqRecord(read));
        			last2Read = read.getReadName();
        			output2Count += 1;
        		}        		
        	}
        }
                
        output1.close();
        output2.close();
        
        if (output1Count != output2Count) {
        	throw new IllegalStateException("Non-symmetrical read counts found for " + inputSam);
        }
	}
	
	private FastqRecord samReadToFastqRecord(SAMRecord read) {
		byte[] bases = read.getReadBases();
		byte[] qualities = read.getBaseQualities();
		
		if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		FastqRecord fastq = new FastqRecord(read.getReadName(), new String(bases), new String(qualities));
		
		return fastq;
	}
	
	private boolean isFirstInPair(SAMRecord read) {
		boolean isFirstInPair;
		
		if (shouldIdentifyEndByReadId) {
			isFirstInPair = read.getReadName().endsWith("/1");
			
		} else {
			isFirstInPair = read.getFirstOfPairFlag();
		}
		
		return isFirstInPair;
	}
	
		
	public void setShouldIdentifyEndByReadId(boolean shouldIdentifyEndByReadId) {
		this.shouldIdentifyEndByReadId = shouldIdentifyEndByReadId;
	}

	public static void run(String[] args) throws IOException {
		Sam2FastqOptions options = new Sam2FastqOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			Sam2Fastq sam2Fastq = new Sam2Fastq();
			if (options.isPairedEnd()) {
				sam2Fastq.setShouldIdentifyEndByReadId(options.shouldIdEndByReadName());
				sam2Fastq.convert(options.getInputFile(), options.getFastq1(), options.getFastq2());
			} else {
				throw new IllegalArgumentException ("Single end not yet supported for sam2fastq");
			}
		}
	}
}
