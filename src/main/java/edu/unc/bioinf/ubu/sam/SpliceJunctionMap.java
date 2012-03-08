package edu.unc.bioinf.ubu.sam;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads a file containing splice junctions creating a map of splice junction 
 * identifiers keyed by splice junction objects.
 * 
 * Each line in the splice junction file is represents an identifier in a 
 * format taken from the TCGA GAF:
 * chr1:12227:+,chr1:12595:+
 * 
 * Chromosome and strand should always be the same in the start and stop
 * positions.
 * 
 * Coordinates defined in the file are exclusive
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SpliceJunctionMap {
    
    private Map<SpliceJunction, String> junctions = new HashMap<SpliceJunction, String>();
    private List<SpliceJunction> junctionList = new ArrayList<SpliceJunction>();

    public SpliceJunctionMap(String junctionFile) throws FileNotFoundException , IOException {
        init(junctionFile);
    }
    
    public List<SpliceJunction> getAllSpliceJunctions() {
//        List<SpliceJunction> junctionList = new ArrayList<SpliceJunction>();
//        junctionList.addAll(junctions.keySet());
//        Collections.sort(junctionList);
        return junctionList;
    }
    
    public void init(String junctionFile) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(junctionFile));
        
        String line = reader.readLine();
        
        while (line != null) {
            String[] range = line.split(",");
            
            if (range.length != 2) {
                throw new IllegalArgumentException("Invalid splice junction format: [" + line + "]");
            }
            
            String[] start = range[0].split(":");
            String[] stop = range[1].split(":");
            
            if ((start.length != 3) || (stop.length != 3)) {
                throw new IllegalArgumentException("Invalid splice junction format: [" + line + "]");
            }
            
            String startChromosome = start[0];
            int startPos           = toInt(start[1], line);
            String startStrand     = start[2];
            
            String stopChromosome = stop[0];
            int stopPos           = toInt(stop[1], line);
            String stopStrand     = stop[2];
            
            if (!(startChromosome.equals(stopChromosome)) || 
                !(startStrand.equals(stopStrand))) {
                throw new IllegalArgumentException("Invalid splice junction format: [" + line + "]");
            }
            
            if (startPos > stopPos) {
                int temp = startPos;
                startPos = stopPos;
                stopPos = temp;
            }
            
            // Adjust coordinates to make them inclusive
            // TODO - Double check negative strand coordinates
            SpliceJunction junction = new SpliceJunction(startChromosome, startPos+1, stopPos-1);
            junctions.put(junction, line);
            junctionList.add(junction);
//            junctions.put(new SpliceJunction(startChromosome, startPos, stopPos), line);
                        
            line = reader.readLine();
        }
        
        reader.close();
    }
    
    private int toInt(String value, String line) {
        try {
            return Integer.valueOf(value);
        } catch (NumberFormatException e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Invalid splice junction format: [" + line + "]");
        }
    }
    
    public String getJunctionKey(SpliceJunction junction) {
        return junctions.get(junction);
    }
    
    public static void main(String[] args) throws Exception {
        SpliceJunctionMap map = new SpliceJunctionMap("/home/lisle/gaf/splice_junctions.txt");
        
        System.out.println("Size: " + map.junctions.size());
    }
}
