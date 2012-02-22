package edu.unc.bioinf.ubu.sam;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections.MultiHashMap;
import org.apache.commons.collections.MultiMap;

/**
 * This class can be used to read a bed file and cache Isoforms information.
 * Isoforms can be retrieved by chromosome and genome position.
 * 
 * TODO(lmose): This class could use better naming throughout.
 * @author lmose
 */
public class BedReader {
    
    private int readOffset;
    
    private Map<String, MultiHashMap> idxMap = new HashMap<String, MultiHashMap>();
    
    private Map<String, Isoform> isoformMap = new HashMap<String, Isoform>();
    
    private List<Coordinate> getExonCoordinates(int isoformStart, String exonLengths, String exonOffsets) {
        List<Coordinate> coordinates = new ArrayList<Coordinate>();
        
        String[] lengths = exonLengths.split(",");
        String[] offsets = exonOffsets.split(",");
        
        if (lengths.length != offsets.length) {
            throw new IllegalArgumentException("Invalid exon length offset pairing: [" + exonLengths + "] [" + exonOffsets + "]");
        }
        
        for (int i = 0; i<lengths.length; i++) {
            // Exon coordinates cached relative to Genome.
            int start = isoformStart;
            start += Integer.parseInt(offsets[i]);
            int stop = start + Integer.parseInt(lengths[i]) - 1;
            
            coordinates.add(new Coordinate(start, stop));
        }
        
        return coordinates;
    }

    /**
     * Build index for fast lookup of isoforms matching a genomic location. 
     */
    public void buildReadToIsoformIndex(String filename, int readOffset) throws Exception {
        
        this.readOffset = readOffset;
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        
        String line = reader.readLine();
        
        int cnt = 1;
        
        while (line != null) {
            String[] fields = line.split("\\s+");
            String chromosome = fields[0];
            // Bed is 0 based, SAM is 1 based - Add one to cache 1 based coordinates. 
            int start = Integer.valueOf(fields[1]) + 1;
            int end = Integer.valueOf(fields[2]) + 1;
            String isoformId = fields[3];
            String strand = fields[5];
            String exonLengths = fields[10];
            String exonOffsets = fields[11];
            
            int startIdx = start - (start % readOffset);
            int endIdx   = end + (readOffset - (end % readOffset));
            
            for (int idx = startIdx; idx <= endIdx; idx += readOffset) {
                getCoordinateMap(chromosome).put(idx, isoformId);
            }
            
            if ((cnt++ % 10000) == 0) {
                //break;
                System.out.println("Loaded " + cnt + " isoforms.");
            }
            
            Isoform isoform = new Isoform(isoformId, new Coordinate(start, end), strand, getExonCoordinates(start, exonLengths, exonOffsets));
            isoformMap.put(isoformId, isoform);
            
            line = reader.readLine();
        }
    }
    
    /**
     * Returns all currently cached Isoforms.
     */
    public Collection<Isoform> getAllIsoforms() {
        return Collections.unmodifiableCollection(isoformMap.values());
    }
    
    private MultiMap getCoordinateMap(String chromosome) {
        MultiMap map = idxMap.get(chromosome);
        
        if (map == null) {
            map = new MultiHashMap();
            idxMap.put(chromosome, (MultiHashMap) map);
        }
        
        return map;
    }
    
    private Collection<String> getPotentialIsoformsExact(String chromosome, int genomeStart) {
        Collection<String> isoforms = Collections.EMPTY_LIST;
        
        MultiMap map = idxMap.get(chromosome);
        if (map != null) {
            Collection<String> isoformCol = (Collection<String>) map.get(genomeStart);
            if (isoformCol != null) {
                isoforms = isoformCol;
            }
        }
        
        return isoforms;
    }
    
    /**
     * Returns a Collection of isoforms that potentially match the input chromosome and genome position
     */
    public Collection<String> getPotentialIsoforms(String chromosome, int genomePosition) {
        Set<String> isoforms = new HashSet<String>();
        
        isoforms.addAll(getPotentialIsoformsExact(chromosome, genomePosition - (genomePosition % readOffset)));
        isoforms.addAll(getPotentialIsoformsExact(chromosome, genomePosition + (readOffset - (genomePosition % readOffset))));
        
        return Collections.unmodifiableCollection(isoforms);
    }
    
    /**
     * Returns the Isoform matching the input isoformId
     */
    public Isoform getIsoform(String isoformId) {
        return isoformMap.get(isoformId);
    }
        
    // Just used for testing.  Users of this class should instantiate and use directly.
    public static void main(String[] args) throws Exception {

        long s = System.currentTimeMillis();
        
        BedReader rdr = new BedReader();
        
        rdr.buildReadToIsoformIndex("/home/lisle/data/coord_convert/ucsc_known_gene_bed.txt", 200);
        
        long e = System.currentTimeMillis();
        
        System.out.println("elapsed secs: " + (e-s)/1000);
        
        System.out.println(rdr.getPotentialIsoforms("chr22", 30000000));
    }
}
