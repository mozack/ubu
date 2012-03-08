package edu.unc.bioinf.ubu.sam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Counts total and unique hits for reads against a given gene.
 * Input BAM is expected to be aligned to transcriptome.
 * 
 * TODO: Either merge with GeneReadCounter or discard.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class GeneReadCounterBWA {
	
	private Map<String, Long> totalGeneCounts = new HashMap<String, Long>(); 
	private Map<String, Long> uniqueGeneCounts = new HashMap<String, Long>();
	
	private IsoformGeneMap isoformGeneMap = new IsoformGeneMap();
	
	private static final Long ONE = (long) 1;
	private static final Long ZERO = (long) 0;

	public void count(IsoformGeneMap isoformGeneMap, String inputSam, String outputFile) throws IOException {
		this.isoformGeneMap = isoformGeneMap;
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
		
		countGenes(inputSam);
		outputCounts(writer);
		
		writer.close();
	}
	
	private void countGenes(String inputFile) {
		
        File file = new File(inputFile);
        
        SAMFileReader reader = new SAMFileReader(file);
        reader.setValidationStringency(ValidationStringency.SILENT);

        int cnt = 0;
        for (SAMRecord read : reader) {
        	if (!read.getReadUnmappedFlag()) {
        		Set<String> genes = new HashSet<String>();
        		
        		String isoform = read.getReferenceName();
        		String gene = isoformGeneMap.getGene(isoform);
        		genes.add(gene);
        		
        		int x0 = 1;
        		Integer x0Tag = read.getIntegerAttribute("X0");
        		if (x0Tag != null) {
        			x0 = x0Tag;
        		}
        		
        		String xaTag = read.getStringAttribute("XA");

        		if (xaTag != null) {
	        		String[] isos = xaTag.split(";");
	        		for (String iso : isos) {
	        			String[] fields = iso.split(",");
	        			if (fields.length > 0) {
	        				String isoId = fields[0];
	        				String geneId = isoformGeneMap.getGene(isoId);
	        				genes.add(geneId);
	        			}
	        		}
        		}
        		
        		if ((genes.size() == 1) && (x0 < 100)) {
        			incrementCount(uniqueGeneCounts, gene);
        		}
        		
        		incrementCount(totalGeneCounts, gene);
        	}
        	
            if ((cnt++ % 1000000) == 0) {
            	System.out.println("Processed " + cnt + " reads.");
            }
        }
	}
	
	private void outputCounts(BufferedWriter writer) throws IOException {
		List<String> genes = isoformGeneMap.getSortedGeneList();
		
		for (String gene : genes) {
			StringBuffer line = new StringBuffer();
			
			long totalCount = getCount(totalGeneCounts, gene);
			
			long uniqueCount = getCount(uniqueGeneCounts, gene);
			
			line.append(gene);
			line.append('\t');
			line.append(totalCount);
			line.append('\t');
			line.append(uniqueCount);
			
//			line.append('\t');
//			line.append(totalCount - uniqueCount);
			
			line.append('\n');
			
			writer.write(line.toString());
		}
	}
	
	private void incrementCount(Map<String, Long> counts, String gene) {
		Long count = counts.get(gene);
		if (count == null) {
			counts.put(gene, ONE);
		} else {
			counts.put(gene, count + 1);
		}
	}
	
	private long getCount(Map<String, Long> counts, String gene) {
		Long count = counts.get(gene);
		if (count == null) {
			count = ZERO;
		}
		return count;
	}
	
	public static void main(String[] args) throws IOException {
				
		String isoformGeneFile = args[0];
		String input = args[1];
		String output = args[2];
		
//		String isoformGeneFile = "/home/lisle/gaf/ref/gaf.knownToLocus";
//		String input = "/home/lisle/data/gene_counts/bwa2.bam";
//		String output = "/home/lisle/data/gene_counts/bwa2.tsv";

		long s = System.currentTimeMillis();
				
		IsoformGeneMap isoformGeneMap = new IsoformGeneMap();
		isoformGeneMap.init(isoformGeneFile);
		
		new GeneReadCounterBWA().count(isoformGeneMap, input, output);
		
		long e = System.currentTimeMillis();
		
		System.out.println("Done.  Elapsed secs: " + (e-s)/1000);
	}
}
