package edu.unc.bioinf.ubu.sam;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Basic map of isoform id to gene id
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class IsoformGeneMap {
	
	public Map<String, String> isoformGeneMap = new HashMap<String, String>();

	public void init(String isoformGeneFile) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(isoformGeneFile));
        
        String line = reader.readLine();
        
        while (line != null) {
        	String[] fields = line.split("\t");
        	if (fields.length == 2) {
        		String gene = fields[0];
        		String isoform = fields[1];
        		
        		if (!isoform.equals(gene)) {
        			isoformGeneMap.put(isoform, gene);
        		}
        	}
            
        	line = reader.readLine();
        	
        }
        
        reader.close();

        System.out.println("Loaded " + isoformGeneMap.keySet().size() + " isoforms");
	}
	
	public List<String> getSortedGeneList() {
		Set<String> geneSet = new HashSet<String>();
		geneSet.addAll(isoformGeneMap.values());
		
		List<String> genes = new ArrayList<String>();
		genes.addAll(geneSet);
		Collections.sort(genes);
		return genes;
	}
	
	public String getGene(String isoform) {
		String gene = isoformGeneMap.get(isoform);
		
		if (gene == null) {
			gene = isoform;
		}
		
		return gene;
	}
}
