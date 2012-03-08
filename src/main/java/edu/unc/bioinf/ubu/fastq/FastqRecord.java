package edu.unc.bioinf.ubu.fastq;

import java.util.Arrays;

/**
 * Representation of a single Fastq record.
 * Some code here may be specific to paired end processing.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqRecord {
    
    public static final int NUM_LINES = 4;
    
    private String[] lines = new String[NUM_LINES];
    
    public FastqRecord(String[] lines) {
        if (lines.length != NUM_LINES) {
            throw new IllegalArgumentException("Invalid number of lines for FastqRecord: [" + Arrays.toString(lines) + "]");
        }
        this.lines = lines;
    }
    
    public String getId() {
        return lines[0];
    }
    
    public String[] getLines() {
        return lines;
    }
    
    /**
     * Returns the portion of the id string leading up to "/"
     */
    public String getBaseId() {
        int slashIdx = getId().indexOf("/");
        int spaceIdx = getId().indexOf(" ");
        
        if ((slashIdx == -1) && (spaceIdx == -1)) {
            return getId();
        }
        
        int idx = -1;
        if (slashIdx == -1) {
            idx = spaceIdx;
        } else if (spaceIdx == -1) {
            idx = slashIdx;
        } else {
            idx = spaceIdx < slashIdx ? spaceIdx : slashIdx;
        }
        
        return getId().substring(0, idx);
    }
    
    /**
     * Returns true if this FastqRecord has the same base id as the input FastqRecord 
     */
    public boolean hasSameBaseId(FastqRecord rec) {
        return rec != null && this.getBaseId().equals(rec.getBaseId());
    }
    
    public String toString() {
        return Arrays.toString(lines);
    }
    
    public int hashcode() {
        return Arrays.hashCode(lines);
    }
    
    public boolean equals(Object obj) {
        FastqRecord that = (FastqRecord) obj;
        return Arrays.equals(this.lines, that.lines);
    }
    
    public void stripNonReadInfoInId() {
        int idx = lines[0].indexOf(" ");
        if (idx > 0) {
            lines[0] = lines[0].substring(0, idx);
        }
    }
    
    public void appendToId(String suffix) {
        if (!lines[0].endsWith(suffix)) {
            lines[0] = lines[0] + suffix;
        }
    }
}
