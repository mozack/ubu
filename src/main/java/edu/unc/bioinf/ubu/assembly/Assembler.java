package edu.unc.bioinf.ubu.assembly;

import java.io.BufferedWriter;
import java.io.File;
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

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

import edu.unc.bioinf.ubu.fastq.FastqInputFile;
import edu.unc.bioinf.ubu.fastq.FastqRecord;
import edu.unc.bioinf.ubu.sam.ReadBlock;

public class Assembler {
	
//	static final int KMER_SIZE = 55;
//	static final int MIN_EDGE_FREQUENCY = 10;
//	static final int MIN_NODE_FREQUENCY = 1;
	
	private int kmerSize = 33;
//	private int kmerSize = 77;
//	private int minEdgeFrequency = 25;
//	private int minNodeFrequncy = 25;
	
	private int minEdgeFrequency = 15;
	private int minNodeFrequncy = 15;

//	private int minContigLength = 101;
	private int minContigLength = 101;
	//private double minEdgeRatio = .51;
	private double minEdgeRatio = .05;
	
	private int minMergeSize = 25;
	
	private FastqInputFile fastq = new FastqInputFile();
	
	private Map<String, Node> nodes = new HashMap<String, Node>();
	
	private Set<Node> rootNodes = new HashSet<Node>();
	
	private List<Contig> contigs = new ArrayList<Contig>();
	
	private BufferedWriter writer;
	
	private SAMFileHeader samHeader;
	
	//private Map<String, SAMRecord> updatedReads = new HashMap<String, SAMRecord>();
	Set<SAMRecord> updatedReads = new HashSet<SAMRecord>();
	
	//private Aligner aligner = new Aligner("/home/lisle/reference/chr5/chr5.fa");
	
	private Aligner aligner = new Aligner("/home/lisle/reference/chr17/chr17.fa");
	
	//private Aligner aligner = new Aligner("/home/lisle/reference/chrX/chrX.fa");
	
	public void assemble(String inputSam, String outputPrefix) throws Exception {
		String contigsFasta = outputPrefix + "_contigs.fasta";
		String contigsSam = outputPrefix + "_contigs.sam";
		String readsBam = outputPrefix + "_reads.bam";
		
		System.out.println("Assembling contigs");
		assembleContigs(inputSam, contigsFasta);
		
		System.out.println("Aligning contigs");
		alignContigs(contigsFasta, contigsSam);
		
		System.out.println("Adjusting reads");
		adjustReads(contigsSam);
		
		System.out.println("Writing adjusted reads");
		outputReads(readsBam);
		
		System.out.println("Done.");
	}
	
	public void assembleContigs(String inputSam, String output) throws FileNotFoundException, IOException {
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        samHeader = reader.getFileHeader();

		writer = new BufferedWriter(new FileWriter(output, false));
		
		int numRecs = 0;
		
		for (SAMRecord read : reader) {
			addToGraph(read);
			numRecs++;
		}
		
		System.out.println("Num records: " + numRecs);
		System.out.println("Num nodes: " + nodes.size());
		
		printEdgeCounts();
		
		filterLowFrequencyEdges();
		filterLowFrequencyNodes();
		
		identifyRootNodes();
		
		buildContigs();
//		mergeContigs();
		outputContigs();
		
		writer.close();
		reader.close();
	}
	
	private void adjustReads(String contigSam) {
		Map<String, Contig> contigMap = new HashMap<String, Contig>();
		for (Contig contig : contigs) {
			contigMap.put(contig.getDescriptor(), contig);
		}
		
        SAMFileReader reader = new SAMFileReader(new File(contigSam));
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        for (SAMRecord contigRead : reader) {
        	
        	List<ReadBlock> contigReadBlocks = ReadBlock.getReadBlocks(contigRead);
        	Contig contig = contigMap.get(contigRead.getReadName());
        	List<ReadPosition> readPositions = contig.getFilteredReadPositions();
        	for (ReadPosition readPosition : readPositions) {
        		//TODO: Won't handle multi-mappers
        		
//        		int position = contigRead.getAlignmentStart() + readPosition.getPosition();
        		
//        		if (position != readPosition.getRead().getAlignmentStart()) {
//        			System.out.println("Different position: " + readPosition.getRead().getReadName());
//        		}
        		
//    			SAMRecord read = readPosition.getRead();
//    			read.setAlignmentStart(position);
    			SAMRecord updatedRead = updateReadAlignment(contigReadBlocks, readPosition);
    			if (updatedRead != null) {
    				updatedReads.add(updatedRead);
    			}
    			//updatedReads.put(read.getReadName(), read);
        		
    			/*
        		SAMRecord updatedRead = updatedReads.get(readPosition.getRead().getReadName());
        		if (updatedRead == null) {
        			SAMRecord read = readPosition.getRead();
        			read.setAlignmentStart(position);
        			updatedReads.put(read.getReadName(), read);
        		} else if (updatedRead.getAlignmentStart() != position) {
        			System.out.println("Skipping multiple position for read: " + updatedRead.getReadName());
        		}
        		*/
        	}
        }
	}
	
    private SAMRecord cloneRead(SAMRecord read) {
        try {
            return (SAMRecord) read.clone();
        } catch (CloneNotSupportedException e) {
            // Infamous "this should never happen" comment here.
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
	
	SAMRecord updateReadAlignment(List<ReadBlock> contigReadBlocks, ReadPosition orig) {
		List<ReadBlock> blocks = new ArrayList<ReadBlock>();
		SAMRecord read = cloneRead(orig.getRead());
		
		int contigPosition = orig.getPosition();
		int accumulatedLength = 0;
		
		// read block positions are one based
		// ReadPosition is zero based
		
		for (ReadBlock contigBlock : contigReadBlocks) {
			if ((contigBlock.getReadStart() + contigBlock.getReferenceLength()) >= orig.getPosition() + 1) {
				ReadBlock block = contigBlock.getSubBlock(accumulatedLength, contigPosition, read.getReadLength() - accumulatedLength);
				
				//TODO: Investigate how this could happen
				if (block.getLength() != 0) {
					blocks.add(block);
					
					if (block.getType() != CigarOperator.D) {
						accumulatedLength += block.getLength();
					}
					
					if (accumulatedLength > read.getReadLength()) {
						throw new IllegalStateException("Accumulated Length: " + accumulatedLength + " is greater than read length: " + read.getReadLength());
					}
					
					if (accumulatedLength == read.getReadLength()) {
						break;
					}
				}
			}
		}
		
		//TODO: Investigate how this could happen.
		if (blocks.size() > 0) {
			int newAlignmentStart = blocks.get(0).getReferenceStart();
			String newCigar = ReadBlock.toCigarString(blocks);
			
			read.setCigarString(newCigar);
			read.setAlignmentStart(newAlignmentStart);
		} else {
			read = null;
		}
		
		return read;
	}
	
	private void outputReads(String readsBam) {
		
		samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(samHeader,
                true, new File(readsBam));
 
        for (SAMRecord read : updatedReads) {
        	out.addAlignment(read);
        }
        
        out.close();
	}
	
	private void alignContigs(String contigFile, String contigSam) throws InterruptedException, IOException {
		aligner.align(contigFile, contigSam);
	}
	
	private void _assemble(String inputFastq, String output) throws FileNotFoundException, IOException {
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
//		mergeContigs();
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
			contig.setDescriptor("contig" + count++ + "_" + contig.getDescriptor());
			writer.append(">" + contig.getDescriptor() + "\n");
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
	
	private void buildContigs() {
		System.out.println("Num starting nodes: " + rootNodes.size());
		
		for (Node node : rootNodes) {
			//StringBuffer contig = new StringBuffer();
			Contig contig = new Contig();
			Set<Node> visitedNodes = new HashSet<Node>();
			Counts counts = new Counts();
			buildContig(node, visitedNodes, contig, counts);
		}
	}
	
	private void processContigTerminus(Node node, Counts counts, Contig contig) {
		
		if (!counts.isTerminatedAtRepeat()) {
			// We've reached the terminus, append the remainder of the node.
			contig.append(node, node.getSequence());
		}
		
		if (contig.getSequence().length() >= minContigLength) {
			contig.setDescriptor(counts.toString());
			contigs.add(contig);
		}
	}
	
	private void buildContig(Node node, Set<Node> visitedNodes, Contig contig, Counts counts) {
		if (visitedNodes.contains(node)) {
			counts.setTerminatedAtRepeat(true);
			processContigTerminus(node, counts, contig);
		} else {
			visitedNodes.add(node);
			
			Collection<Edge> edges = node.getToEdges();
			
			if (edges.isEmpty()) {
				processContigTerminus(node, counts, contig);

			} else {
				// Append current character
				contig.append(node, Character.toString(node.getSequence().charAt(0)));
				
				// Create a new contig branch for each edge
				for (Edge edge : edges) {
					counts.incrementEdgeCounts(edge.getCount());
					Contig contigBranch = new Contig(contig);
					Set<Node> visitedNodesBranch = new HashSet<Node>(visitedNodes);
					buildContig(edge.getTo(), visitedNodesBranch, contigBranch, (Counts) counts.clone());
				}				
			}			
		}
	} 
	
	// Merge contigs that overlap with < kmerSize bases
	// This addresses "smallish" gaps in the graph
	private void mergeContigs() {
		
		if (minMergeSize > kmerSize) {
			List<Contig> updatedContigs = new ArrayList<Contig>(contigs);
			
			int mergedCount = 0;
			
			for (Contig contig1 : contigs) {
				boolean isMerged = false;
				
				for (Contig contig2 : updatedContigs) {
					if (contig1 != contig2) {
						int overlapIdx = getOverlapIndex(contig1.getSequence(), contig2.getSequence());
						
						if (overlapIdx > -1) {
							contig2.prependSequence(contig1.getDescriptor(), contig1.getSequence());
							isMerged = true;
						}
					}
				}
				
				if (isMerged) {
					updatedContigs.remove(contig1);
					mergedCount += 1;
				}
			}
			
			this.contigs = updatedContigs;
			
			System.out.println("Merged: " + mergedCount + " overlapping contigs.");
		}
	}
	
	private int getOverlapIndex(String s1, String s2) {
		int strLenDiff = s2.length() - s1.length();
		
		// Default start to 0 or the length of s2 from the end of s1
		int start = strLenDiff > 0 ? strLenDiff : 0;
		
		// If minMergeSize from end of s1 is beyond start, update start
		start = Math.max(start, s1.length()-minMergeSize);
		
		for (int i=start; i<s1.length(); i++) {
			if (s2.startsWith(s1.substring(i))) {
				return i;
			}
		}
		
		return -1;
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
	
	private void addToGraph(SAMRecord read) {
		Node node = addToGraph(read.getReadString());
		
		if (node != null) {
			node.addStartingRead(read);
		}
	}
	
	private Node addToGraph(String sequence) {
		Node prev = null;
		Node firstNode = null;
		
		for (int i=0; i<=sequence.length()-kmerSize; i++) {
			String kmer = sequence.substring(i, i+kmerSize);
			Node node = nodes.get(kmer);
			if (node == null) {
				node = new Node(kmer, sequence);
				nodes.put(kmer, node);
			} else {
				node.incrementCount();
			}
			
			if (prev != null) {
				prev.addToEdge(node);
			}
			
			if (firstNode == null) {
				firstNode = node;
			}
			
			prev = node;
		}
		
		return firstNode;
	}
	
	public static void main(String[] args) throws Exception {
		long s = System.currentTimeMillis();
		
		Assembler ayc = new Assembler();
		
//		ayc.assemble("/home/lisle/ayc/case0/normal_7576572_7577692.fastq", "/home/lisle/ayc/case0/normal_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/case0/tumor_7576572_7577692.fastq", "/home/lisle/ayc/case0/tumor_33_05.fasta");

		
		//ayc.assemble("/home/lisle/ayc/tp53.fastq", "/home/lisle/ayc/tp53.fasta");
//		ayc.assemble("/home/lisle/ayc/run2/normal_7576572_7577692.fastq", "/home/lisle/ayc/run4/normal_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/run2/tumor_7576572_7577692.fastq", "/home/lisle/ayc/run4/tumor_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/case1/normal.fastq", "/home/lisle/ayc/case1/normal_33_05.fasta");
//		ayc.assemble("/home/lisle/ayc/case1/tumor.fastq", "/home/lisle/ayc/case1/tumor_33_05.fasta");
		
//		ayc.assemble("/home/lisle/ayc/case2/normal.fastq", "/home/lisle/ayc/case2/deeper/rcnormal_19_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/tumor.fastq", "/home/lisle/ayc/case2/deeper/rctumor_19_02.fasta");
		
//		ayc.assemble("/home/lisle/ayc/case2/normal.fastq", "/home/lisle/ayc/case2/deeper/normal_77_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/tumor.fastq", "/home/lisle/ayc/case2/deeper/tumor_77_02.fasta");
		
//		ayc.assemble("/home/lisle/ayc/case2/normal.fastq", "/home/lisle/ayc/case2/deeper/normal_77_02B.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/tumor.fastq", "/home/lisle/ayc/case2/deeper/tumor_77_02B.fasta");

//		ayc.assemble("/home/lisle/ayc/case1/round2/tumor.bam", "/home/lisle/ayc/case1/round2/re_tumor");
		
//		ayc.assemble("/home/lisle/ayc/case2/round2/case2_tumor.bam", "/home/lisle/ayc/case2/round2/ra_tumor");
		
		ayc.assemble("/home/lisle/ayc/case0/round2/case0_tumor.bam", "/home/lisle/ayc/case0/round2/ra_tumor");
		
//		ayc.assemble("/home/lisle/ayc/case2/realigned/ra_normal.fastq", "/home/lisle/ayc/case2/realigned/normal_33_02.fasta");
//		ayc.assemble("/home/lisle/ayc/case2/realigned/ra_tumor.fastq", "/home/lisle/ayc/case2/realigned/tumor_33_02.fasta");


		
		long e = System.currentTimeMillis();
		
		System.out.println("Elapsed secs: " + (e-s)/1000);
	}
}
