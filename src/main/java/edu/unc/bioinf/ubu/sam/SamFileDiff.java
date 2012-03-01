package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Diffs 2 SAM (or BAM) files.
 * 
 * @author lmose
 */
public class SamFileDiff {
	
	private SAMRecord cachedRead1;
	private SAMRecord cachedRead2;
	
	private SAMRecordComparator comparator = new SAMRecordComparator();
	
	private List<SAMRecord> cachedReadList1;
	private List<SAMRecord> cachedReadList2;
	private Iterator<List<SAMRecord>> iter1;
	private Iterator<List<SAMRecord>> iter2;

	public void diff(String samInputFileName1, String samInputFileName2, String samOutputFileName1, String samOutputFileName2) {
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
        	        	
        	int compare = readList1.get(0).getSAMString().compareTo(readList2.get(0).getSAMString());
        	if (compare < 0) {
        		addAlignments(out1, readList1);
        		this.cachedReadList2 = readList2;
        	} else if (compare > 0) {
        		addAlignments(out2, readList2);
        		this.cachedReadList1 = readList1;
        	} else {
        		diffReadLists(readList1, readList2, out1, out2);
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
	
	static class SAMRecordComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord read1, SAMRecord read2) {
			return read1.getSAMString().compareTo(read2.getSAMString());
		}
	}
	
	public static void main(String[] args) {
		String in1  = args[0];
		String in2  = args[1];
		String out1 = args[2];
		String out2 = args[3];

//		String in1  = "/home/lisle/data/sam_diff/small/1.sam";
//		String in2  = "/home/lisle/data/sam_diff/small/2.sam";
//		String out1 = "/home/lisle/data/sam_diff/small/out1.sam";
//		String out2 = "/home/lisle/data/sam_diff/small/out2.sam";

		System.out.println("Input 1: " + in1);
		System.out.println("Input 2: " + in2);
		System.out.println("Output 1: " + out1);
		System.out.println("Output 2: " + out2);
		
		long s = System.currentTimeMillis();
		
		new SamFileDiff().diff(in1, in2, out1, out2);
		
		long e = System.currentTimeMillis();
		
		System.out.println("Done.  Elapsed secs: " + (e-s)/1000);
	}
}
