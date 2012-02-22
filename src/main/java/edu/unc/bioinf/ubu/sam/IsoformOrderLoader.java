package edu.unc.bioinf.ubu.sam;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Reads a fasta file and caches ordering of isoforms in a Map
 * 
 * @author lmose
 */
public class IsoformOrderLoader {

    private Map<String, Integer> isoformOrderMap;
    private static final char ISOFORM_INDICATOR = '>';
    
    public void loadOrdering(String filename) throws FileNotFoundException, IOException {
        isoformOrderMap = new HashMap<String, Integer>();
        
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        
        String line = reader.readLine();
        
        int cnt = 1;
        
        while (line != null) {

            if (line.length() < 2) {
                throw new IllegalArgumentException("Invalid entry [" + line + "] in file [" + filename + "]");
            }
            
            if (line.charAt(0) == ISOFORM_INDICATOR) {
                String isoformId = line.substring(1, line.length());
                isoformOrderMap.put(isoformId, cnt);
                cnt++;
            }
            
            line = reader.readLine();
        }
        
        System.out.println("Ordering loaded for " + cnt + " isosforms.");
        
        reader.close();
    }
    
    public Integer getOrder(String isoformId) {
        if (isoformOrderMap == null) {
            throw new IllegalStateException("Isoform order map not initialized.");
        }
        
        Integer order = isoformOrderMap.get(isoformId);
        
        if (order == null) {
            order = Integer.MAX_VALUE;
        }
        
        return order;
    }
    
    public void clearCache() {
        // This is a somewhat large map, so make it eligible for GC
        isoformOrderMap = null;
    }
}
