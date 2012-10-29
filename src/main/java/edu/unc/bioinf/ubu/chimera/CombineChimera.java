package edu.unc.bioinf.ubu.chimera;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.unc.bioinf.ubu.sam.SamMultiMappingReader;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class CombineChimera {

	/**
	 * Roll chimeric reads up into a single indel if they are on the same chromosome and strand
	 * Input BAM must be sorted by name.
	 */
	public void combine(String input, String output, int minMappingQuality) {
		SamMultiMappingReader reader = new SamMultiMappingReader(input);
		
		SAMFileHeader header = reader.getFileHeader();
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriter outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, true, new File(output));

		for (List<SAMRecord> readList : reader) {
			
			boolean isCombined = false;
			
			List<SAMRecord> prunedReadList = pruneLowMappingQualityReads(readList, minMappingQuality);
			
			if (prunedReadList.size() == 2) {
				SAMRecord read1 = prunedReadList.get(0);
				SAMRecord read2 = prunedReadList.get(1);
				
				// Look for same reference, strand and multipart Cigar
				if ((read1.getReferenceName().equals(read2.getReferenceName())) &&
					(read1.getReadNegativeStrandFlag() == read2.getReadNegativeStrandFlag()) &&
					(read1.getCigarLength() >= 2) &&
					(read2.getCigarLength() >= 2)) {
					
					SAMRecord combinedRead = combineChimericReads(read1, read2);
					if (combinedRead != null) {
						isCombined = true;
						
	//					outputReadsBam.addAlignment(combinedRead);
						outputRead(combinedRead, isCombined, outputReadsBam);
					}
				}
			}
			
			if (!isCombined) {
				for (SAMRecord origRead : readList) {
//					outputReadsBam.addAlignment(origRead);
					outputRead(origRead, isCombined, outputReadsBam);
				}
			}
		}
		
		outputReadsBam.close();
		reader.close();
	}
	
	private void outputRead(SAMRecord read, boolean isCombined, SAMFileWriter out) {
		try {
			out.addAlignment(read);
		} catch (NullPointerException e) {
			System.out.println("isCombined: " + isCombined);
			System.out.println(read);
			System.out.println(read.getReadName());
			System.out.println(read.getSAMString());
			throw e;
		}
	}
	
	private SAMRecord combineChimericReads(SAMRecord read1, SAMRecord read2) {
		SAMRecord combinedRead = null;
		
		Cigar left = null;
		Cigar right = null;
		int leftPos = 0;
		int rightPos = 0;
		
		// Look for S at end of left side and S at start of right side
		// Left side position should be less that right side position
		if ((read1.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) &&
			(read2.getCigar().getCigarElement(read2.getCigarLength()-1).getOperator() == CigarOperator.S) &&
			(read2.getAlignmentStart() <  read1.getAlignmentStart())) {
			right = read1.getCigar();
			rightPos = read1.getAlignmentStart();
			left = read2.getCigar();
			leftPos = read2.getAlignmentStart();
		} else if ((read2.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) &&
				   (read1.getCigar().getCigarElement(read1.getCigarLength()-1).getOperator() == CigarOperator.S) &&
				   (read1.getAlignmentStart() < read2.getAlignmentStart())) {
			right = read2.getCigar();
			rightPos = read2.getAlignmentStart();
			left = read1.getCigar();			
			leftPos = read1.getAlignmentStart();
		}
		
		if ((left != null) && (right != null)) {
			
			List<CigarElement> leftElements = new ArrayList<CigarElement>();
			List<CigarElement> rightElements = new ArrayList<CigarElement>();
			
			// Drop trailing S on left side
			for (int i=0; i<left.numCigarElements()-1; i++) {
				leftElements.add(left.getCigarElement(i));
			}

			// Drop leading S on right side
			for (int i=1; i<right.numCigarElements(); i++) {
				rightElements.add(right.getCigarElement(i));
			}
			
			int totalLength = this.getTotalLength(leftElements, rightElements);
			
			int trimmedElemLength = 100;
			
			// If total element length is longer than the read, then trim first element
			// on the right side of the indel (this is likely a deletion??)
			if (totalLength > read1.getReadLength()) {
				// Trim from the left side of the leftmost block on the right of the indel
				int trimLength = totalLength - read1.getReadLength();
				trimmedElemLength = trimLeftmostElement(rightElements, trimLength);
				rightPos += trimLength; 
			}
			
			int leftLength = this.getTotalLength(leftElements);
			// left end is exclusive
			int leftEnd = leftPos + leftLength;
			
			// If left side starts with S, the leftPos is set to the start of the next element.
			// Adjust the left end accordingly
			if (leftElements.get(0).getOperator() == CigarOperator.S) {
				leftEnd -= leftElements.get(0).getLength();
			}
			
			// If the end of the left side of the read is greater than the start of the right,
			// Trim the first element on the right side of the read (this is likely an insertion??).
			if (leftEnd > rightPos) {
				int rightShift = leftEnd - rightPos;
				trimmedElemLength = trimLeftmostElement(rightElements, rightShift);
				rightPos = leftEnd;
			}
			
			// Create representation of indel element
			CigarOperator indelType;
			int indelLength;
			if (rightPos > leftEnd) {
				indelType = CigarOperator.D;
				indelLength = rightPos - (leftEnd);
			} else {
				indelType = CigarOperator.I;
				indelLength = read1.getReadLength() - getTotalLength(leftElements, rightElements);
			}
			
			if (trimmedElemLength > 0) {
				// Re-create element list with left side elements, indel, and right side elements
				List<CigarElement> elements = new ArrayList<CigarElement>();
				elements.addAll(leftElements);
				elements.add(new CigarElement(indelLength, indelType));
				elements.addAll(rightElements);
				
				// Create combined read
				combinedRead = cloneRead(read1);
				combinedRead.setAlignmentStart(leftPos);
				combinedRead.setCigar(new Cigar(elements));
				combinedRead.setMappingQuality((read1.getMappingQuality() + read2.getMappingQuality()) / 2);
			}
		}

		return combinedRead;
	}
	
	private SAMRecord cloneRead(SAMRecord read) {
		try {
			return (SAMRecord) read.clone();
		} catch (CloneNotSupportedException e) {
			// Infamous "this should never happen" comment here.
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	private int trimLeftmostElement(List<CigarElement> rightElements, int trimLength) {
		CigarElement toTrim = rightElements.get(0);
		int newLength = toTrim.getLength() - trimLength;
		CigarElement replacement = new CigarElement(newLength, toTrim.getOperator());
		rightElements.set(0, replacement);
		
		return newLength;
	}
	
	private int getTotalLength(List<CigarElement> elements) {
		int total = 0;
		
		for (CigarElement element : elements) {
			if (element.getOperator() != CigarOperator.D) {
				total += element.getLength();
			}
		}
		
		return total;
	}
	
	private int getTotalLength(List<CigarElement> left, List<CigarElement> right) {
		return getTotalLength(left) + getTotalLength(right);
	}
	
	private List<SAMRecord> pruneLowMappingQualityReads(List<SAMRecord> reads, int minMappingQuality) {
		
		List<SAMRecord> pruned = new ArrayList<SAMRecord>();
		
		for (SAMRecord read : reads) {
			if (read.getMappingQuality() >= minMappingQuality) {
				pruned.add(read);
			}
		}
		
		return pruned;
	}
	
	public static void main(String[] args) {
//		String in = "/home/lmose/dev/ayc/long_indels/32I.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/32I_out.sam";
		
//		String in = "/home/lmose/dev/ayc/long_indels/25I.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/25I_out.sam";
		
//		String in = "/home/lmose/dev/ayc/long_indels/25D.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/25D_out.sam";

//		String in = "/home/lmose/dev/ayc/long_indels/100D.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/100D_out.sam";
//		int mapq = 2;

//		String in = "/home/lmose/dev/ayc/long_indels/100D_sc.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/100D_sc_out.sam";
//		int mapq = 2;

//		String in = "/home/lmose/dev/ayc/long_indels/s411.sam";
//		String out = "/home/lmose/dev/ayc/long_indels/s411_out.sam";
//		int mapq = 2;

		
		String in = args[0];
		String out = args[1];
		int mapq = Integer.parseInt(args[2]);
		
		CombineChimera cc = new CombineChimera();
		cc.combine(in, out, mapq);
	}
}
