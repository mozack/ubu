package edu.unc.bioinf.ubu.fastq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Handles reading FastqRecords from file.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqInputFile {
    
    private static final int CACHING_DISABLED = -1;
    
    private BufferedReader reader;
    private Map<Integer, FastqRecord> records = new HashMap<Integer, FastqRecord>();
    private int recordNum = 0;
    private int maxCachedLines;
    
    public void init(String filename, int maxCachedLines) throws FileNotFoundException {
        openFile(filename);
        this.maxCachedLines = maxCachedLines;
    }
    
    public void init (String filename) throws FileNotFoundException {
        init(filename, CACHING_DISABLED);
    }
    
    private void openFile(String filename) throws FileNotFoundException {
        reader = new BufferedReader(new FileReader(filename));
    }

    public FastqRecord getNextRecord() throws IOException {
        String[] lines = new String[FastqRecord.NUM_LINES];
        
        for (int i=0; i<4; i++) {
            String line = reader.readLine();
            if (line == null) {
                return null;
            }
            lines[i] = line;
        }
        
        return new FastqRecord(lines);
    }
    
    /**
     *  Retrieves the record 
     */
    public FastqRecord getRecord(int idx) throws IOException {
        
        if (isCachingDisabled()) {
            throw new UnsupportedOperationException("Retrieving record by index not supported when caching is disabled.");
        }
        
        if (idx < 1) {
            throw new IllegalArgumentException("Read index cannot be less than 1.  Index: [" + idx + "]");
        }
        
        FastqRecord record = null;
        
        // Attempt to get record from cache
        if (records.containsKey(idx)) {
            return records.get(idx);
        }
        
        // Can't "read back" before the cache
        if (idx <= recordNum) {
            throw new UnsupportedOperationException("Cannot read back to an uncached record. Index: [" + idx + "]");
        }
        
        // Read forward to the desired record and update the cache
        while (idx > recordNum) {
            record = getNextRecord();
            recordNum++;
            // Cache the newly read record and purge stale records.
            // i.e. Put record 101 and remove record 1.  
            records.put(recordNum, record);
            records.remove(recordNum - maxCachedLines);
        }
        
        // The cache should never be bigger than maxCachedLines
        assert(records.size() <= maxCachedLines);
        
        return record;
    }
    
    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }
    
    private boolean isCachingDisabled() {
        return maxCachedLines == CACHING_DISABLED;
    }
}
