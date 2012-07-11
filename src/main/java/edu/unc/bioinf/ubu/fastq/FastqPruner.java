package edu.unc.bioinf.ubu.fastq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class FastqPruner {
	
	private Set<String> readsToFilter = new HashSet<String>();

	public void prune(String inputSam, String inFastq1, String inFastq2, String outFastq1, String outFastq2, int bit) 
		throws FileNotFoundException, IOException {
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        System.out.println("Loading reads to filter.");
        
        for (SAMRecord read : reader) {
        	if ((read.getFlags() & bit) == bit) {
        		readsToFilter.add(getBaseId(read.getReadName()));
        	}
        }
        
        System.out.println("Filter size: " + readsToFilter.size());
        
        reader.close();
        
        System.out.println("Pruning " + inFastq1);
        int count1 = prune(inFastq1, outFastq1);
        
        System.out.println("Pruning " + inFastq2);
        int count2 = prune(inFastq2, outFastq2);
        
        System.out.println("Num records written to " + outFastq1 + " : " + count1);
        System.out.println("Num records written to " + outFastq2 + " : " + count2);
        
        System.out.println("Done.");
	}
	
	private int prune(String inFastq, String outFastq) throws FileNotFoundException, IOException {
        FastqInputFile in = new FastqInputFile();
        in.init(inFastq);
        
        FastqOutputFile out = new FastqOutputFile();
        out.init(outFastq);

        int count = 0;
        
        FastqRecord rec = in.getNextRecord();
        
        while (rec != null) {
        	if (!readsToFilter.contains((rec.getBaseId()))) {
        		out.write(rec);
        		count += 1;
        	}
        	
        	rec = in.getNextRecord();
        }
        
        in.close();
        out.close();
        
        return count;
	}
	
	// TODO: dupe of FastqRecord.getBaseId
    private String getBaseId(String id) {
        int slashIdx = id.indexOf("/");
        int spaceIdx = id.indexOf(" ");
        
        if ((slashIdx == -1) && (spaceIdx == -1)) {
            return id;
        }
        
        int idx = -1;
        if (slashIdx == -1) {
            idx = spaceIdx;
        } else if (spaceIdx == -1) {
            idx = slashIdx;
        } else {
            idx = spaceIdx < slashIdx ? spaceIdx : slashIdx;
        }
        
        return id.substring(0, idx);
    }
    
    public static void main(String[] args) throws Exception {
    	String inputSam  = args[0];
    	String inFastq1  = args[1];
    	String inFastq2  = args[2];
    	String outFastq1 = args[3];
    	String outFastq2 = args[4];
    	int bit          = Integer.parseInt(args[5]);
    	
    	long s = System.currentTimeMillis();
    	
    	new FastqPruner().prune(inputSam, inFastq1, inFastq2, outFastq1, outFastq2, bit);
    	
    	long e = System.currentTimeMillis();
    	
    	System.out.println("Elapsed: " + (e-s)/1000);
    }
}
