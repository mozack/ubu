package edu.unc.bioinf.ubu.fastq;

import java.io.File;
import java.io.IOException;
import java.util.List;

import edu.unc.bioinf.ubu.sam.ReadBlock;
import edu.unc.bioinf.ubu.sam.ReverseComplementor;

import net.sf.samtools.CigarOperator;
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
	private String end1Suffix;
	private String end2Suffix;
	private boolean shouldHandleMapspliceFusions;
	
	private static final String MAPSLICE_FUSION_TAG = "ZF";
	private static final String MAPSLICE_FUSION     = "FUS";

	/**
	 * Convert the input paired end SAM/BAM file into 2 fastq files.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 * Note: muti-mappers that do not contain /1 or /2 in the read name may not be handled properly due to sort order.
	 */
	public void convert(String inputSam, String outputFastq1, String outputFastq2) throws IOException {
		String last1Read = "";
		String last2Read = "";
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        output1 = new FastqOutputFile();
        output1.init(outputFastq1);
        output2 = new FastqOutputFile();
        output2.init(outputFastq2);
        
        int output1Count = 0;
        int output2Count = 0;
        int lineCnt = 0;
        
        SAMRecord fusionRead = null;
        
        for (SAMRecord read : reader) {
        	if (isFirstInPair(read)) {
        		if (!read.getReadName().equals(last1Read)) {
        			last1Read = read.getReadName();
        			if (!isFusion(read)) {
	        			output1.write(samReadToFastqRecord(read));
	        			output1Count += 1;
	        			fusionRead = null;
        			} else {
        				fusionRead = read;
        			}
        		} else {
        			if ((isFusion(read)) && (fusionRead != null) && (read.getAttribute(MAPSLICE_FUSION_TAG).equals(fusionRead.getAttribute(MAPSLICE_FUSION_TAG)))) {
        				output1.write(fusionToFastqRecord(fusionRead, read));
        				output1Count += 1;
        				fusionRead = null;
        			}
        		}
        	} else if (isSecondInPair(read)) {
        		if (!read.getReadName().equals(last2Read)) {
        			last2Read = read.getReadName();
        			if (!isFusion(read)) {
	        			output2.write(samReadToFastqRecord(read));
	        			output2Count += 1;
	        			fusionRead = null;
        			} else {
        				fusionRead = read;
        			}
        		} else {
        			if ((isFusion(read)) && (fusionRead != null) && (read.getAttribute(MAPSLICE_FUSION_TAG).equals(fusionRead.getAttribute(MAPSLICE_FUSION_TAG)))) {
        				output2.write(fusionToFastqRecord(fusionRead, read));
        				output2Count += 1;
        				fusionRead = null;
        			}        			
        		}
        	}
        	
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.out.println("record: " + lineCnt);
            }
        }
                
        output1.close();
        output2.close();
        
        if (output1Count != output2Count) {
        	throw new IllegalStateException("Non-symmetrical read counts found for " + inputSam + ".  Your reads may not be paired properly.");
        }
	}
	
	private boolean isFusion(SAMRecord read) {
		boolean isFusion = false;
		
		if (shouldHandleMapspliceFusions) {
			String fusionTag = (String) read.getAttribute(MAPSLICE_FUSION_TAG);
			isFusion = (fusionTag != null) && (fusionTag.contains(MAPSLICE_FUSION));
		}
		
		return isFusion;
	}
	
	/**
	 * Convert the input SAM/BAM file into a single fastq file.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFastq) throws IOException {
		String last1Read = "";
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        output1 = new FastqOutputFile();
        output1.init(outputFastq);
        int lineCnt = 0;
        
        for (SAMRecord read : reader) {
    		if (!read.getReadName().equals(last1Read)) {
    			output1.write(samReadToFastqRecord(read));
    			last1Read = read.getReadName();
    		}
    		
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.out.println("record: " + lineCnt);
            }
        }
                
        output1.close();
	}
	
	private FastqRecord samReadToFastqRecord(SAMRecord read) {
		String bases = read.getReadString();
		String qualities = read.getBaseQualityString();
		
		if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		FastqRecord fastq = new FastqRecord("@" + read.getReadName(), bases, qualities);
		
		return fastq;
	}
	
	private FastqRecord fusionToFastqRecord(SAMRecord read1, SAMRecord read2) {
		
		SAMRecord doner = null;
		SAMRecord accepter = null;
		int donerLength = -1;
		int accepterStart = -1;
		
		List<ReadBlock> read1Blocks = ReadBlock.getReadBlocks(read1);
		List<ReadBlock> read2Blocks = ReadBlock.getReadBlocks(read2);
		
		ReadBlock read1FirstBlock = read1Blocks.get(0);
		ReadBlock read1LastBlock  = read1Blocks.get(read1Blocks.size()-1);
		ReadBlock read2FirstBlock = read2Blocks.get(0);
		ReadBlock read2LastBlock  = read2Blocks.get(read2Blocks.size()-1);
		
		if ((read1FirstBlock.getType() == CigarOperator.S) &&
			(read2LastBlock.getType() == CigarOperator.S)) {
			
			accepter = read1;
			accepterStart = read1FirstBlock.getLength() + 1;

			doner = read2;
			donerLength = read2.getReadLength() - read2LastBlock.getLength();
			
		} else if ((read1LastBlock.getType() == CigarOperator.S) &&
				   (read2FirstBlock.getType() == CigarOperator.S)) {
			
			doner = read1;
			donerLength = read1.getReadLength() - read1LastBlock.getLength();
			
			accepter = read2;
			accepterStart = read2FirstBlock.getLength();
		}
		
		String donerBases = doner.getReadString().substring(0, donerLength);
		if (doner.getReadNegativeStrandFlag()) {
			donerBases = reverseComplementor.reverseComplement(donerBases);
		}
		
		String accepterBases = accepter.getReadString().substring(accepterStart);
		if (accepter.getReadNegativeStrandFlag()) {
			accepterBases = reverseComplementor.reverseComplement(accepterBases);
		}
		
		String bases = donerBases + accepterBases;
		
		// Assumption: Qualities are in the original format
		String qualities = doner.getBaseQualityString();
		
		// Read name should be same for both reads.
		String readName = doner.getReadName();
		
		FastqRecord fastq = new FastqRecord("@" + readName, bases, qualities);
		
		return fastq;
	}
	
	private boolean isFirstInPair(SAMRecord read) {
		boolean isFirstInPair;
		
		if (shouldIdentifyEndByReadId) {
			isFirstInPair = read.getReadName().endsWith(end1Suffix);
			
		} else {
			isFirstInPair = read.getFirstOfPairFlag();
		}
		
		return isFirstInPair;
	}
	
	private boolean isSecondInPair(SAMRecord read) {
		boolean isSecondInPair;
		
		if (shouldIdentifyEndByReadId) {
			isSecondInPair = read.getReadName().endsWith(end2Suffix);
			
		} else {
			isSecondInPair = read.getSecondOfPairFlag();
		}
		
		return isSecondInPair;
	}
	
	public void setEndSuffixes(String end1Suffix, String end2Suffix) {
		this.shouldIdentifyEndByReadId = true;
		this.end1Suffix = end1Suffix;
		this.end2Suffix = end2Suffix;
	}
	
	public void setShouldHandleMapspliceFusions(boolean shouldHandleMapspliceFusions) {
		this.shouldHandleMapspliceFusions = shouldHandleMapspliceFusions;
	}
	
	public static void run(String[] args) throws IOException {
		Sam2FastqOptions options = new Sam2FastqOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			long s = System.currentTimeMillis();
			System.out.println("sam2fastq starting");
			
			Sam2Fastq sam2Fastq = new Sam2Fastq();
			sam2Fastq.setShouldHandleMapspliceFusions(options.shouldHandleMapspliceFusions());
			
			if (options.isPairedEnd()) {
				
				if (options.shouldIdEndByReadName()) {
					sam2Fastq.setEndSuffixes(options.getEnd1Suffix(), options.getEnd2Suffix());
				}
				
				sam2Fastq.convert(options.getInputFile(), options.getFastq1(), options.getFastq2());
			} else {
				//throw new IllegalArgumentException("Single end is not currently supported");
				System.out.println("WARNING: Single end not thoroughly tested!\n");
				sam2Fastq.convert(options.getInputFile(), options.getFastq1());
			}
			
			long e = System.currentTimeMillis();
			System.out.println("sam2fastq done.  Elapsed secs: " + (e-s)/1000);
		}
	}
	
	public static void main(String[] args) throws Exception {
		String[] argz =
			"--in /home/lisle/sam2fastq/fusion.sam --fastq1 /home/lisle/sam2fastq/fus1.fastq --fastq2 /home/lisle/sam2fastq/fus2.fastq --end1 /1 --end2 /2 --mfusion".split(" ");
		
		run(argz);
	}
}
