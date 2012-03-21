package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang.StringUtils;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Diffs 2 SAM (or BAM) files.
 * The input files must be sorted by read
 * Input file sorting should be done via samtools.
 * (OS level sort may produce different results).
 * <p>
 * 2 output SAM (or BAM) files are produced containing reads contained
 * in one file but not the other.  The entire SAM string is used to identify
 * a read.
 * <p>
 * Handles multi-mapped reads (as produced by Mapsplice).
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamFileDiff {
	
	private SAMRecord cachedRead1;
	private SAMRecord cachedRead2;
	
	private SAMRecordComparator comparator = new SAMRecordComparator();
	
	private List<SAMRecord> cachedReadList1;
	private List<SAMRecord> cachedReadList2;
	private Iterator<List<SAMRecord>> iter1;
	private Iterator<List<SAMRecord>> iter2;
	
	private boolean isReadIdComparisonOnly;

	public void diff(String samInputFileName1, String samInputFileName2, String samOutputFileName1, String samOutputFileName2) {
		
		System.out.println("in1:\t" + samInputFileName1);
		System.out.println("in2:\t" + samInputFileName2);
		System.out.println("out1:\t" + samOutputFileName1);
		System.out.println("out2:\t" + samOutputFileName2);
		
		SamMultiMappingReader in1 = new SamMultiMappingReader(samInputFileName1);
		SamMultiMappingReader in2 = new SamMultiMappingReader(samInputFileName2);
		
        final SAMFileWriter out1 = new SAMFileWriterFactory().makeSAMOrBAMWriter(in1.getFileHeader(),
                true, new File(samOutputFileName1));
        
        final SAMFileWriter out2 = new SAMFileWriterFactory().makeSAMOrBAMWriter(in2.getFileHeader(),
                true, new File(samOutputFileName2));

        iter1 = in1.iterator();
        iter2 = in2.iterator();
        
        while ((hasNextList1()) && (hasNextList2())) {
        	List<SAMRecord> readList1 = getNextList1();
        	List<SAMRecord> readList2 = getNextList2();
        	        	
        	int compare = compareReadNames(readList1.get(0), readList2.get(0));
        	
        	if (compare < 0) {
        		addAlignments(out1, readList1);
        		this.cachedReadList2 = readList2;
        	} else if (compare > 0) {
        		addAlignments(out2, readList2);
        		this.cachedReadList1 = readList1;
        	} else {
        		if (!isReadIdComparisonOnly) {
        			diffReadLists(readList1, readList2, out1, out2);
        		}
        	}
        }

        while (hasNextList1()) {
        	addAlignments(out1, getNextList1());
        }
        
        while (hasNextList2()) {
        	addAlignments(out2, getNextList2());
        }
        
        out1.close();
        out2.close();
	}
	
	public void setReadIdComparisonOnly(boolean isReadIdComparisonOnly) {
		this.isReadIdComparisonOnly = isReadIdComparisonOnly;
	}
	
	private void addAlignments(SAMFileWriter out, List<SAMRecord> reads) {
		for (SAMRecord read : reads) {
			out.addAlignment(read);
		}
	}
	
	private boolean hasNextList1() {
		return (cachedReadList1 != null) || iter1.hasNext();
	}
	
	private boolean hasNextList2() {
		return (cachedReadList2 != null) || iter2.hasNext();
	}
	
	private List<SAMRecord> getNextList1() {
		List<SAMRecord> list;
		
		if (cachedReadList1 != null) {
			list = cachedReadList1;
			cachedReadList1 = null;
		} else {
			list = iter1.next();
		}
		
		return list;
	}
	
	private List<SAMRecord> getNextList2() {
		List<SAMRecord> list;
		
		if (cachedReadList2 != null) {
			list = cachedReadList2;
			cachedReadList2 = null;
		} else {
			list = iter2.next();
		}
		
		return list;
	}
	
	private boolean hasNextRead1(Iterator<SAMRecord> iter) {
		return this.cachedRead1 != null || iter.hasNext();
	}
	
	private boolean hasNextRead2(Iterator<SAMRecord> iter) {
		return this.cachedRead2 != null || iter.hasNext();
	}
	
	private SAMRecord getNextRead1(Iterator<SAMRecord> iter) {
		SAMRecord read = null;
		
		if (this.cachedRead1 != null) {
			read = cachedRead1;
			cachedRead1 = null;
		} else {
			read = iter.next();
		}
		
		return read;
	}
	
	private SAMRecord getNextRead2(Iterator<SAMRecord> iter) {
		SAMRecord read = null;
		
		if (this.cachedRead2 != null) {
			read = cachedRead2;
			cachedRead2 = null;
		} else {
			read = iter.next();
		}
		
		return read;
	}

	
	private void diffReadLists(List<SAMRecord> readList1, List<SAMRecord> readList2, SAMFileWriter out1, SAMFileWriter out2) {
		
		Collections.sort(readList1, comparator);
		Collections.sort(readList2, comparator);
				
		Iterator<SAMRecord> iter1 = readList1.iterator();
		Iterator<SAMRecord> iter2 = readList2.iterator();
		
		// Iterate through reads in parallel
		//while ((iter1.hasNext()) && (iter2.hasNext())) {
        while ((hasNextRead1(iter1)) && (hasNextRead2(iter2))) {
        	SAMRecord read1 = getNextRead1(iter1);
        	SAMRecord read2 = getNextRead2(iter2);
        	
        	int compare = read1.getSAMString().compareTo(read2.getSAMString());
        	
        	if (compare < 0) {
        		// There is an extra read in the first bam file.
        		out1.addAlignment(read1);
        		this.cachedRead2 = read2;
        	} else if (compare > 0) {
        		// There is an extra read in the second bam file.
        		out2.addAlignment(read2);
        		this.cachedRead1 = read1;
        	}
        }
        
        // Add leftovers at end (if any)
        while (hasNextRead1(iter1)) {
        	SAMRecord read1 = getNextRead1(iter1);
        	out1.addAlignment(read1);
        }
        
        while (hasNextRead2(iter2)) {
        	SAMRecord read2 = getNextRead2(iter2);
        	out2.addAlignment(read2);
        }
        
        this.cachedRead1 = null;
        this.cachedRead2 = null;
	}
	
	private int compareReadNames(SAMRecord read1, SAMRecord read2) {
		int compare = 0;
		
		String[] name1 = read1.getReadName().split("[:/]");
		String[] name2 = read2.getReadName().split("[:/]");
		
		int idx = 0;
		while ((idx < name1.length) && (idx < name2.length) && (compare == 0)) {
			if (StringUtils.isNumeric(name1[idx]) && StringUtils.isNumeric(name2[idx])) {
				int field1 = Integer.parseInt(name1[idx]);
				int field2 = Integer.parseInt(name2[idx]);
				
				compare = field1 - field2;
			} else {
				compare = name1[idx].compareTo(name2[idx]);
			}
			
			idx++;
		}
		
		return compare;
	}
	
	static class SAMRecordComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord read1, SAMRecord read2) {
			return read1.getSAMString().compareTo(read2.getSAMString());
		}
	}
	
	public static void run(String[] args) {
		
		SamFileDiffOptions options = new SamFileDiffOptions();
		options.parseOptions(args);
		if (options.isValid()) {
			
			long s = System.currentTimeMillis();
			
			SamFileDiff diff = new SamFileDiff();
			
			diff.setReadIdComparisonOnly(options.isReadIdComparisonOnly());
			
			diff.diff(options.getInput1File(), options.getInput2File(),
					options.getOutput1File(), options.getOutput2File());
			
			long e = System.currentTimeMillis();
			
			System.out.println("Done.  Elapsed secs: " + (e-s)/1000);
		}
	}
}
