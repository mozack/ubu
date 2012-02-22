package edu.unc.bioinf.ubu.sam;

import net.sf.samtools.SAMRecord;

/**
 * A Pair of SAMRecord's
 * 
 * @author lmose
 */
public class ReadPair {
    private SAMRecord read1;
    private SAMRecord read2;
    
    ReadPair(SAMRecord read1, SAMRecord read2) {
        this.read1 = read1;
        this.read2 = read2;
    }
    
    public SAMRecord getRead1() {
        return read1;
    }
    
    public SAMRecord getRead2() {
        return read2;
    }
    
    public String toString() {
        String r1 = read1 != null ? read1.getReadName() : "null";
        String r2 = read2 != null ? read2.getReadName() : "null";
        return "read1: " + r1 + ", read2: " + r2;
    }
}
