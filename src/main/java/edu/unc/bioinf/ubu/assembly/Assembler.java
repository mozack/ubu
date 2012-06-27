package edu.unc.bioinf.ubu.assembly;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.unc.bioinf.ubu.fastq.FastqInputFile;
import edu.unc.bioinf.ubu.fastq.FastqRecord;
import edu.unc.bioinf.ubu.sam.ReverseComplementor;

public class Assembler {
	
//	static final int KMER_SIZE = 55;
//	static final int MIN_EDGE_FREQUENCY = 10;
//	static final int MIN_NODE_FREQUENCY = 1;
	
	private int kmerSize = 77;
//	private int minEdgeFrequency = 25;
//	private int minNodeFrequncy = 25;
	
	private int minEdgeFrequency = 10;
	private int minNodeFrequncy = 10;

//	private int minContigLength = 101;
	private int minContigLength = 101;
	//private double minEdgeRatio = .51;
	private double minEdgeRatio = .02;
	
	private FastqInputFile fastq = new FastqInputFile();
	
	private Map<String, Node> nodes = new HashMap<String, Node>();
	
	private Set<Node> rootNodes = new HashSet<Node>();
	
	private List<Contig> contigs = new ArrayList<Contig>();
	
	private BufferedWriter writer; 
	
	public void assemble(String inputFastq, String output) throws FileNotFoundException, IOException {
		fastq.init(inputFastq);
		writer = new BufferedWriter(new FileWriter(output, false));
		
		FastqRecord rec = fastq.getNextRecord();
		int numRecs = 0;
		
		while (rec != null) {
			String sequence = rec.getSequence();
//			addToGraph(new ReverseComplementor().reverseComplement(sequence));
			addToGraph(sequence);
			
			rec = fastq.getNextRecord();
			numRecs++;
		}
		
		System.out.println("Num records: " + numRecs);
		System.out.println("Num nodes: " + nodes.size());
		
		printEdgeCounts();
		
		filterLowFrequencyEdges();
		filterLowFrequencyNodes();
		
		identifyRootNodes();
		
		buildContigs();
		outputContigs();
		
		writer.close();
	}
	
	public void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
	}
	
	public void setMinContigLength(int minContigLength) {
		this.minContigLength = minContigLength;
	}

	public void setMinEdgeFrequency(int minEdgeFrequency) {
		this.minEdgeFrequency = minEdgeFrequency;
	}

	public void setMinNodeFrequncy(int minNodeFrequncy) {
		this.minNodeFrequncy = minNodeFrequncy;
	}
	
	public void setMinEdgeRatio(double minEdgeRatio) {
		this.minEdgeRatio = minEdgeRatio;
	}

	private void filterLowFrequencyNodes() {
		List<Node> nodesToFilter = new ArrayList<Node>();
		
		for (Node node : nodes.values()) {
			if (node.getCount() < minNodeFrequncy) {
				nodesToFilter.add(node);
			}
		}
		
		List<Edge> edgesToFilter = new ArrayList<Edge>();
		
		for (Node node : nodesToFilter) {
			edgesToFilter.addAll(node.getToEdges());
			edgesToFilter.addAll(node.getFromEdges());
		}
		
		for (Edge edge : edgesToFilter) {
			edge.remove();
		}
		
		for (Node node : nodesToFilter) {
			nodes.remove(node.getSequence());
		}
	}
	
	private void filterLowFrequencyEdges() {
		
		Set<Edge> edgesToFilter = new HashSet<Edge>();
		
		for (Node node : nodes.values()) {
			for (Edge edge : node.getToEdges()) {
				if (edge.getCount() < minEdgeFrequency) {
					edgesToFilter.add(edge);
				}
			}
			
			edgesToFilter.addAll(node.getInfrequentEdges(minEdgeRatio));
		}
		
		for (Edge edge : edgesToFilter) {
			edge.remove();
		}
		

	}
	
	private void outputContigs() throws IOException {
		
		int count = 0;
		
		for (Contig contig : contigs) {
			writer.append(">contig" + count++ + "_" + contig.getDescriptor() + "\n");
			writer.append(contig.getSequence());
			writer.append("\n");
		}
	}
	
	private void identifyRootNodes() {
		for (Node node : nodes.values()) {
			if (node.isRootNode()) {
				rootNodes.add(node);
			}
		}
	}
	
	private boolean containsMultipleSources(Set<Node> nodes) {
		String prevSource = null;
		
		for (Node node : nodes) {
			if ((prevSource != null) && (!node.getSource().equals(prevSource))) {
				return true;
			}
			
			prevSource = node.getSource();
		}
		
		return false;
	}
	
	private void buildContigs() {
		System.out.println("Num starting nodes: " + rootNodes.size());
		
		for (Node node : rootNodes) {
			StringBuffer contig = new StringBuffer();
			Set<Node> visitedNodes = new HashSet<Node>();
			Counts counts = new Counts();
			buildContig(node, visitedNodes, contig, counts);
			
//			if (containsMultipleSources(visitedNodes)) {
		}
	}
	
	private void processContigTerminus(Node node, Counts counts, StringBuffer contig) {
		
		if (!counts.isTerminatedAtRepeat()) {
			// We've reached the terminus, append the remainder of the node.
			contig.append(node.getSequence());
		}
		
		String contigStr = contig.toString();
		
		if (contigStr.length() >= minContigLength) {
			contigs.add(new Contig(counts.toString(), contig.toString()));
		}
	}
	
	private void buildContig(Node node, Set<Node> visitedNodes, StringBuffer contig, Counts counts) {
		if (visitedNodes.contains(node)) {
			counts.setTerminatedAtRepeat(true);
			processContigTerminus(node, counts, contig);
		} else {
			visitedNodes.add(node);
			
			//TODO: Handle heterozygous case
//			Edge edge = node.getMostCommonEdge();
			
			//This should be uneccessary now with pre-filtering
			//List<Edge> edges = node.getFrequentToEdges(minEdgeRatio);
			Collection<Edge> edges = node.getToEdges();
			
			if (edges.isEmpty()) {
				processContigTerminus(node, counts, contig);

			} else {
				// Append current character
				contig.append(node.getSequence().charAt(0));
				
				// Create a new contig branch for each edge
				for (Edge edge : edges) {
					counts.incrementEdgeCounts(edge.getCount());
					StringBuffer contigBranch = new StringBuffer(contig);
					Set<Node> visitedNodesBranch = new HashSet<Node>(visitedNodes);
					buildContig(edge.getTo(), visitedNodesBranch, contigBranch, (Counts) counts.clone());
				}				
			}			
		}
	}
	
	
	private void printEdgeCounts() {
		long[] edgeCounts = new long[nodes.size()];
		List<Integer> edgeSizes = new ArrayList<Integer>();
		
		int idx = 0;
		for (Node node : nodes.values()) {
			edgeCounts[idx++] = node.getToEdges().size();
			
			for (Edge edge : node.getToEdges()) {
				edgeSizes.add(edge.getCount());
			}
		}
		
		Arrays.sort(edgeCounts);
		System.out.println("Median edge count: " + edgeCounts[edgeCounts.length/2]);
		System.out.println("Max edge count: " + edgeCounts[edgeCounts.length-1]);
		System.out.println("Min edge count: " + edgeCounts[0]);
		Integer[] sizes = edgeSizes.toArray(new Integer[edgeSizes.size()]);
		Arrays.sort(sizes);
		System.out.println("Median edge size: " + sizes[sizes.length/2]);
		System.out.println("Max edge size: " + sizes[sizes.length-1]);
		System.out.println("Min edge size: " + sizes[0]);
	}
	
	private void addToGraph(String sequence) {
		Node prev = null;
		
		for (int i=0; i<=sequence.length()-kmerSize; i++) {
			String kmer = sequence.substring(i, i+kmerSize);
			Node node = nodes.get(kmer);
			if (node == null) {
				node = new Node(kmer, sequence);
				nodes.put(kmer, node);
//				rootNodes.add(node);
			} else {
				node.incrementCount();
			}
			
			if (prev != null) {
//				rootNodes.remove(node);
				prev.addToEdge(node);
			}
			
			prev = node;
		}
	}
	
	public static void main(String[] args) throws Exception {
		long s = System.currentTimeMillis();
		
		Assembler ayc = new Assembler();
		//ayc.assemble("/home/lisle/ayc/tp53.fastq", "/home/lisle/ayc/tp53.fasta");
//		ayc.assemble("/home/lisle/ayc/run2/normal_7576572_7577692.fastq", "/home/lisle/ayc/run4/normal_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/run2/tumor_7576572_7577692.fastq", "/home/lisle/ayc/run4/tumor_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/case1/normal.fastq", "/home/lisle/ayc/case1/normal_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/case1/tumor.fastq", "/home/lisle/ayc/case1/tumor_33_05.fasta");
		
//		ayc.assemble("/home/lisle/ayc/case2/normal.fastq", "/home/lisle/ayc/case2/deeper/rcnormal_19_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/tumor.fastq", "/home/lisle/ayc/case2/deeper/rctumor_19_02.fasta");
		
		ayc.assemble("/home/lisle/ayc/case2/normal.fastq", "/home/lisle/ayc/case2/deeper/normal_77_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/tumor.fastq", "/home/lisle/ayc/case2/deeper/tumor_77_02.fasta");

		
//		ayc.assemble("/home/lisle/ayc/case2/realigned/ra_normal.fastq", "/home/lisle/ayc/case2/realigned/normal_33_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/realigned/ra_tumor.fastq", "/home/lisle/ayc/case2/realigned/tumor_33_02.fasta");


		
		long e = System.currentTimeMillis();
		
		System.out.println("Elapsed secs: " + (e-s)/1000);
	}
}
