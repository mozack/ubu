package edu.unc.bioinf.ubu.assembly;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import edu.unc.bioinf.ubu.fastq.Sam2Fastq;
import edu.unc.bioinf.ubu.gtf.Feature;
import edu.unc.bioinf.ubu.gtf.GtfLoader;
import edu.unc.bioinf.ubu.sam.ReadBlock;
import edu.unc.bioinf.ubu.sam.ReverseComplementor;

public class ReAligner {

	private SAMFileHeader samHeader;

	private SAMFileWriter outputReadsBam;

	private long startMillis;

	private List<Feature> regions;

	private String regionsGtf;

	private String tempDir;

	private String reference;
	
	private String referenceDir;
	
	private int minContigMapq;

	private AssemblerSettings assemblerSettings;
	
	private int numThreads;
	
	private int allowedMismatchesFromContig = 0;
	
	private boolean shouldReprocessUnaligned = true;
	
	private List<ReAlignerRunnable> threads = new ArrayList<ReAlignerRunnable>();
	
	private List<SAMRecord> unalignedReads = new ArrayList<SAMRecord>();
	
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	
	public void reAlign(String inputSam, String outputSam) throws Exception {

		System.out.println("input: " + inputSam);
		System.out.println("output: " + outputSam);
		System.out.println("regions: " + regionsGtf);
		System.out.println("reference: " + reference);
		System.out.println("reference: " + referenceDir);
		System.out.println("working dir: " + tempDir);
		System.out.println("num threads: " + numThreads);
		System.out.println("allowed mismatches: " + allowedMismatchesFromContig);
		System.out.println(assemblerSettings.getDescription());

		startMillis = System.currentTimeMillis();

		init();

		log("Loading target regions");
		loadRegions();

		log("Reading Input SAM Header");
		getSamHeader(inputSam);

		log("Initializing output SAM File");
		initOutputFile(outputSam);
		
		log("Iterating over regions");
		for (Feature region : regions) {
			//processRegion(region, inputSam);
			log("Spawning thread for: " + region.getDescriptor());
			spawnRegionThread(region, inputSam);
		}
		
		log("Waiting for all threads to complete");
		waitForAllThreadsToComplete();
		
		log("Combining contigs");
		String contigFasta = tempDir + "/" + "all_contigs.fasta";
		combineContigs(contigFasta);
		
		processContigs(contigFasta, tempDir, inputSam);
		
		if (shouldReprocessUnaligned) {		
			log("Assembling unaligned reads");
			String unalignedSam = tempDir + "/" + "unaligned_to_contig.bam";
			
			String unalignedDir = tempDir + "/unaligned";
			String unalignedContigFasta = unalignedDir + "/unaligned_contigs.fasta";
			Assembler assem = newAssembler();
			boolean hasContigs = assem.assembleContigs(unalignedSam, unalignedContigFasta, "unaligned");
			
//			String finalUnaligned = unalignedDir + "/" + "unaligned_to_contig.bam";
			
			if (hasContigs) {
				processContigs(unalignedContigFasta, unalignedDir, unalignedSam);
			}
		}
		
		outputReadsBam.close();
		
		/*
		log("Aligning contigs");
		Aligner aligner = new Aligner(reference, numThreads);
		String contigsSam = tempDir + "/" + "all_contigs.sam";
		aligner.align(contigFasta, contigsSam);
		
		log("Cleaning contigs");
		String cleanContigsFasta = tempDir + "/" + "clean_contigs.fasta";
		cleanAndOutputContigs(contigsSam, cleanContigsFasta);
		
		log("Aligning original reads to contigs");
		String alignedToContigSam = tempDir + "/" + "align_to_contig.sam";
		this.alignToContigs(inputSam, alignedToContigSam, cleanContigsFasta);
		
		log("Convert aligned to contig to bam, and sorting bams");
		String alignedToContigBam = tempDir + "/" + "align_to_contig.bam";
		String sortedAlignedToContig = tempDir + "/" + "sorted_aligned_to_contig";
		String sortedOriginalReads = tempDir + "/" + "sorted_original_reads"; 
		runCommand("samtools view -bS " + alignedToContigSam + " -o " + alignedToContigBam);
		runCommand("samtools sort -n " + alignedToContigBam + " " + sortedAlignedToContig);
		runCommand("samtools sort -n " + inputSam + " " + sortedOriginalReads);
		sortedAlignedToContig += ".bam";
		sortedOriginalReads += ".bam";
		
		log("Adjust reads");
		String unaligned = tempDir + "/" + "unaligned_reads.bam";
		adjustReads(sortedOriginalReads, sortedAlignedToContig, unaligned);
		*/

		System.out.println("Done.");
	}
	
	private void processContigs(String contigFasta, String tempDir, String inputSam) throws InterruptedException, IOException {
		log("Aligning contigs");
		Aligner aligner = new Aligner(reference, numThreads);
		String contigsSam = tempDir + "/" + "all_contigs.sam";
		aligner.align(contigFasta, contigsSam);
		
		log("Cleaning contigs");
		String cleanContigsFasta = tempDir + "/" + "clean_contigs.fasta";
		cleanAndOutputContigs(contigsSam, cleanContigsFasta);
		
		log("Aligning original reads to contigs");
		String alignedToContigSam = tempDir + "/" + "align_to_contig.sam";
		this.alignToContigs(inputSam, alignedToContigSam, cleanContigsFasta);
		
		log("Convert aligned to contig to bam, and sorting bams");
		String alignedToContigBam = tempDir + "/" + "align_to_contig.bam";
		String sortedAlignedToContig = tempDir + "/" + "sorted_aligned_to_contig";
		String sortedOriginalReads = tempDir + "/" + "sorted_original_reads"; 
		runCommand("samtools view -bS " + alignedToContigSam + " -o " + alignedToContigBam);
		runCommand("samtools sort -n " + alignedToContigBam + " " + sortedAlignedToContig);
		runCommand("samtools sort -n " + inputSam + " " + sortedOriginalReads);
		sortedAlignedToContig += ".bam";
		sortedOriginalReads += ".bam";
		
		log("Adjust reads");
		String unaligned = tempDir + "/" + "unaligned_to_contig.bam";
		adjustReads(sortedOriginalReads, sortedAlignedToContig, unaligned);
	}
	
	private void combineContigs(String contigFasta) throws IOException, InterruptedException {
		
		for (Feature region : regions) {
			String regionContigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			appendFile(regionContigsFasta, contigFasta);
		}
	}
	
	public synchronized void addThread(ReAlignerRunnable thread) {
		threads.add(thread);
	}
	
	public synchronized void removeThread(ReAlignerRunnable thread) {
		threads.remove(thread);
	}
	
	private synchronized int activeThreads() {
		return threads.size();
	}
	
	private void waitForAvailableThread() throws InterruptedException {
		while (activeThreads() == numThreads) {
			Thread.sleep(500);
		}
	}
	
	private void waitForAllThreadsToComplete() throws InterruptedException {
		long start = System.currentTimeMillis();
		while (activeThreads() > 0) {
			long elapsedSecs = (System.currentTimeMillis() - start) / 1000;
			if ((elapsedSecs % 60) == 0) {
				log("Waiting on " + threads.size() + " threads.");
			}
			Thread.sleep(500);
		}
	}
	
	private void spawnRegionThread(Feature region, String inputSam) throws InterruptedException {
		waitForAvailableThread();
		ReAlignerRunnable thread = new ReAlignerRunnable(this, region, inputSam);
		addThread(thread);
		new Thread(thread).start();
	}
	
	//TODO: Factor out, and use where appropriate
	private void runCommand(String cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.out.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.out.println("cmd time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("cmd exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
	}

	private void runCommand(String[] cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		StringBuffer cmdStr = new StringBuffer();
		for (String substring : cmd) {
			cmdStr.append(substring);
			cmdStr.append(" ");
		}
		
		System.out.println("Running: [" + cmdStr + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.out.println("cmd time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("cmd exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
	}
	
	private void getUnalignedReads(String inputSam) throws InterruptedException, IOException {
		String unalignedBam = tempDir + "/unaligned.bam";
		String unalignedFastq = getUnalignedFastqFile();
		
		String cmd = "samtools view -b -f 0x04 " + inputSam + " -o " + unalignedBam;
		runCommand(cmd);
		
		sam2Fastq(unalignedBam, unalignedFastq);
		
		SAMFileReader reader = new SAMFileReader(new File(unalignedBam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		for (SAMRecord read : reader) {
			if (read.getReadUnmappedFlag()) {
				unalignedReads.add(read);
			}
		}
		reader.close();
	}
	
	private String getUnalignedFastqFile() {
		return tempDir + "/unaligned.fastq";
	}

	public void processRegion(Feature region, String inputSam) throws Exception {
		
		try {
		
			Set<SAMRecord> updatedReads = new HashSet<SAMRecord>();
			
	//		log("Extracting targeted region: " + region.getDescriptor());
			String targetRegionBam = extractTargetRegion(inputSam, region);
			
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			String cleanContigsFasta = tempDir + "/" + region.getDescriptor() + "_clean_contigs.fasta";
			String contigsSam   = tempDir + "/" + region.getDescriptor() + "_contigs.sam";
			String outputBam    = tempDir + "/" + region.getDescriptor() + "_output.bam";
			String targetRegionFastq = tempDir + "/" + region.getDescriptor() + ".fastq";
			String alignedToContigSam = tempDir + "/" + region.getDescriptor() + "_aligned_to_contig.sam";
			String unalignedFastq = getUnalignedFastqFile();
			
			Assembler assem = newAssembler();
			
			assem.assembleContigs(targetRegionBam, contigsFasta, region.getDescriptor());
			
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	private String extractTargetRegion(String inputSam, Feature region)
			throws IOException, InterruptedException {
		String extractFile = tempDir + "/" + region.getDescriptor() + ".bam";

		String location = region.getSeqname() + ":" + region.getStart() + "-"
				+ region.getEnd();

		String cmd = "samtools view -b " + inputSam + " " + location + " -o "
				+ extractFile;

		System.out.println("Running: [" + cmd + "]");

		long s = System.currentTimeMillis();

		Process proc = Runtime.getRuntime().exec(cmd);

		int ret = proc.waitFor();

		long e = System.currentTimeMillis();

		System.out.println("Extract time: " + (e - s) / 1000 + " seconds.");

		if (ret != 0) {
			String stdout = getOutput(proc.getInputStream());
			String stderr = getOutput(proc.getErrorStream());
			System.out.println("Samtools stdout: " + stdout);
			System.out.println("Samtools stderr: " + stderr);

			throw new RuntimeException(
					"Samtools exited with non-zero return code : " + ret);
		}

		return extractFile;
	}

	private String getOutput(InputStream is) throws IOException {
		StringWriter writer = new StringWriter();

		Reader reader = new BufferedReader(new InputStreamReader(is));

		char[] buffer = new char[1024];

		int n;
		while ((n = reader.read(buffer)) != -1) {
			writer.write(buffer, 0, n);
		}

		reader.close();

		return writer.toString();
	}

	private void loadRegions() throws IOException {
		GtfLoader loader = new GtfLoader();
		regions = loader.load(regionsGtf);
	}

	public void setRegionsGtf(String gtfFile) {
		this.regionsGtf = gtfFile;
	}

	private void log(String message) {
		long currMillis = System.currentTimeMillis() - startMillis;
		System.out.println(currMillis + " " + message);
	}

	private void getSamHeader(String inputSam) {
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		reader.setValidationStringency(ValidationStringency.SILENT);

		samHeader = reader.getFileHeader();
		samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

		reader.close();
	}
	
	private void sam2Fastq(String bam, String fastq) throws IOException {
		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.convert(bam, fastq);
	}
	
	private void appendFile(String file1, String file2) throws InterruptedException, IOException {
		String[] cmd = new String[] { "bash", "-c", "cat " + file1 + " >> " + file2 };
		runCommand(cmd);
	}
		
	private void cleanAndOutputContigs(String contigsSam, String cleanContigsFasta) throws IOException {
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(cleanContigsFasta, false));
		
		SAMFileReader contigReader = new SAMFileReader(new File(contigsSam));
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		for (SAMRecord contigRead : contigReader) {
			//TODO: Does this guarantee no alternate alignments?
			if (contigRead.getMappingQuality() > 1) {
				
				removeSoftClips(contigRead);
				
				String bases = contigRead.getReadString();
				
				// Express contigs in forward strand context
				if (contigRead.getReadNegativeStrandFlag()) {
					bases = reverseComplementor.reverseComplement(bases);
				}
				
				//TODO: Safer delimiter?  This assumes no ~ in any read
				contigRead.setReadString("");
				String contigReadStr = contigRead.getSAMString();
				contigReadStr = contigReadStr.replace('\t','~');
				
				String contigName = contigRead.getReadName() + "~" + contigReadStr; 
				
				writer.append(">" + contigName);
				writer.append(bases);
				writer.append("\n");
			}
		}
		contigReader.close();
		
		writer.close();
	}
	
	// Assumes entirely soft clipped reads are filtered prior to here.
	// Checks only the first and last Cigar element
	// Does not adjust qualities
	private void removeSoftClips(SAMRecord read) {
		
		Cigar cigar = read.getCigar();
		
		CigarElement firstElement = cigar.getCigarElement(0);
		CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
		
		if ((firstElement.getOperator() == CigarOperator.S) ||
			(lastElement.getOperator() == CigarOperator.S)) {
		
			Cigar newCigar = new Cigar();
			
			String bases = read.getReadString();
			//String qualities = read.getBaseQualityString();
					
			if (firstElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(firstElement.getLength(), bases.length()-1);
				//qualities = qualities.substring(firstElement.getLength(), qualities.length()-1);
			} else {
				newCigar.add(firstElement);
			}
			
			for (int i=1; i<cigar.numCigarElements()-1; i++) {
				newCigar.add(cigar.getCigarElement(i));
			}
			
			if (lastElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(0, bases.length() - lastElement.getLength() - 1);
				//qualities = qualities.substring(0, qualities.length() - lastElement.getLength() - 1);
			}
			
			read.setCigar(newCigar);
			read.setReadString(bases);
			//read.setBaseQualityString(qualities);
		}
	}
	
	private void alignToContigs(String inputSam, String alignedToContigSam, String contigFasta) throws IOException, InterruptedException {
		// Convert original bam to fastq
		String fastq = tempDir + "/" + "original_reads.fastq";
		sam2Fastq(inputSam, fastq);
		
		// Build contig fasta index
		Aligner contigAligner = new Aligner(contigFasta, numThreads);
		contigAligner.index();
		
		// Align region fastq against assembled contigs
		contigAligner.shortAlign(fastq, alignedToContigSam);
	}
	
	private void adjustReads(String originalReadsSam, String alignedToContigSam, String unalignedSam) throws IOException {
		
		SAMFileWriter unalignedReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(unalignedSam));
		
		SAMFileReader contigReader = new SAMFileReader(new File(alignedToContigSam));
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileReader origReader = new SAMFileReader(new File(originalReadsSam));
		origReader.setValidationStringency(ValidationStringency.SILENT);

		Iterator<SAMRecord> contigIter = contigReader.iterator();
		Iterator<SAMRecord> origIter = origReader.iterator();
		
		SAMRecord cachedContig = null;
		
		SamStringReader samStringReader = new SamStringReader();
		
		while ((contigIter.hasNext() || cachedContig != null) && (origIter.hasNext())) {
			SAMRecord orig = origIter.next();
			SAMRecord read;
			
			if (cachedContig != null) {
				read = cachedContig;
				cachedContig = null;
			} else {
				read = contigIter.next();
			}
			
			if (orig.getReadName().equals(read.getReadName())) {
				//TODO: Check for mismatches.  Smarter CIGAR check.
				if ((read.getCigarString().equals("100M")) && (read.getReadUnmappedFlag() == false)) {
				
					SAMRecord origRead = orig;
					String contigReadStr = read.getReferenceName();
					contigReadStr = contigReadStr.substring(contigReadStr.indexOf('~')+1);
					contigReadStr = contigReadStr.replace('~', '\t');
					SAMRecord contigRead = samStringReader.getRead(contigReadStr);

					List<ReadBlock> contigReadBlocks = ReadBlock.getReadBlocks(contigRead);
					
					ReadPosition readPosition = new ReadPosition(origRead, read.getAlignmentStart()-1, -1);
					SAMRecord updatedRead = updateReadAlignment(contigRead,
							contigReadBlocks, readPosition);
					
					if (updatedRead != null) {
						//TODO: Move into updateReadAlignment ?
						if (updatedRead.getMappingQuality() == 0) {
							updatedRead.setMappingQuality(1);
						}
						
						if (updatedRead.getReadUnmappedFlag()) {
							updatedRead.setReadUnmappedFlag(false);
						}
						
						updatedRead.setReadNegativeStrandFlag(read.getReadNegativeStrandFlag());
						
						// Reverse complement / reverse original read bases and qualities if
						// the original read was unmapped and is now on the reverse strand
						// Originally mapped reads would already be expressed in forward strand context
						// TODO: Do we need to handle forward / reverse strand change.  Is this even possible?
						// TODO: What about reads that align across chromosomes and are tagged as unmapped
						//       Might they already be reverse complemented?
						if ((origRead.getReadUnmappedFlag()) && (read.getReadNegativeStrandFlag())) {
							updatedRead.setReadString(reverseComplementor.reverseComplement(updatedRead.getReadString()));
							updatedRead.setBaseQualityString(reverseComplementor.reverse(updatedRead.getBaseQualityString()));
						}
						
						outputReadsBam.addAlignment(updatedRead);
					}
				}

			} else {
				cachedContig = read;
				if (orig.getReadUnmappedFlag()) {
					unalignedReadsBam.addAlignment(orig);
				}
			}
		}

		unalignedReadsBam.close();
		origReader.close();
		contigReader.close();
		
	}
	
	SAMRecord updateReadAlignment(SAMRecord contigRead,
			List<ReadBlock> contigReadBlocks, ReadPosition orig) {
		List<ReadBlock> blocks = new ArrayList<ReadBlock>();
		SAMRecord read = cloneRead(orig.getRead());

		read.setReferenceName(contigRead.getReferenceName());

		int contigPosition = orig.getPosition();
		int accumulatedLength = 0;

		// read block positions are one based
		// ReadPosition is zero based

		for (ReadBlock contigBlock : contigReadBlocks) {
			if ((contigBlock.getReadStart() + contigBlock.getReferenceLength()) >= orig
					.getPosition() + 1) {
				ReadBlock block = contigBlock.getSubBlock(accumulatedLength,
						contigPosition, read.getReadLength()
								- accumulatedLength);
				
				// If this is an insert, we need to adjust the alignment start
				if ((block.getType() == CigarOperator.I) && (block.getLength() != 0)) {
					contigPosition = contigPosition - (contigBlock.getLength() - block.getLength());
					block.setReferenceStart(block.getReferenceStart() - (contigBlock.getLength() - block.getLength()));
//					block = contigBlock.getSubBlock(accumulatedLength,
//								contigPosition, read.getReadLength()
//								- accumulatedLength);
				}
				
				//TODO: Drop leading and trailing delete blocks

				// TODO: Investigate how this could happen
				if (block.getLength() != 0) {
					blocks.add(block);

					if (block.getType() != CigarOperator.D) {
						accumulatedLength += block.getLength();
					}

					if (accumulatedLength > read.getReadLength()) {
						throw new IllegalStateException("Accumulated Length: "
								+ accumulatedLength
								+ " is greater than read length: "
								+ read.getReadLength());
					}

					if (accumulatedLength == read.getReadLength()) {
						break;
					}
				}
			}
		}

		if (blocks.size() > 0) {
			// If we've aligned past the end of the contig resulting in a short Cigar
			// length, append additonal M to the Cigar
			ReadBlock.fillToLength(blocks, read.getReadLength());
			int newAlignmentStart = blocks.get(0).getReferenceStart();
			String newCigar = ReadBlock.toCigarString(blocks);

			read.setCigarString(newCigar);
			read.setAlignmentStart(newAlignmentStart);
		} else {
			// TODO: Investigate how this could happen.
			read = null;
		}

		return read;
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

	private void initOutputFile(String outputReadsBamFilename) {
		samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

		outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(outputReadsBamFilename));
	}
		
	private Assembler newAssembler() {
		Assembler assem = new Assembler();

		assem.setKmerSize(assemblerSettings.getKmerSize());
		assem.setMinEdgeFrequency(assemblerSettings.getMinEdgeFrequency());
		assem.setMinNodeFrequncy(assemblerSettings.getMinNodeFrequncy());
		assem.setMinContigLength(assemblerSettings.getMinContigLength());
		assem.setMinEdgeRatio(assemblerSettings.getMinEdgeRatio());
		assem.setMaxPotentialContigs(assemblerSettings
				.getMaxPotentialContigs());
		assem.setMinContigRatio(assemblerSettings.getMinContigRatio());
		assem.setMinUniqueReads(assemblerSettings.getMinUniqueReads());

		return assem;
	}

	private void init() {
		File workingDir = new File(tempDir);
		if (workingDir.exists()) {
			if (!workingDir.delete()) {
				throw new IllegalStateException("Unable to delete: " + tempDir);
			}
		}

		if (!workingDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir);
		}
		
		File unalignedTempDir = new File(tempDir + "/unaligned");
		
		if (!unalignedTempDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir + "/unaligned");
		}
	}

	public void setReference(String reference) {
		this.reference = reference;
	}
	
	public void setReferenceDir(String referenceDir) {
		this.referenceDir = referenceDir;
	}

	public void setTempDir(String temp) {
		this.tempDir = temp;
	}

	public void setAssemblerSettings(AssemblerSettings settings) {
		this.assemblerSettings = settings;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public void setMinContigMapq(int minContigMapq) {
		this.minContigMapq = minContigMapq;
	}
	
	public void setAllowedMismatchesFromContig(int allowedMismatchesFromContig) {
		this.allowedMismatchesFromContig = allowedMismatchesFromContig;
	}
	
	public void setShouldReprocessUnaligned(boolean shouldReprocessUnaligned) {
		this.shouldReprocessUnaligned = shouldReprocessUnaligned;
	}

	public static void run(String[] args) throws Exception {
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);

		if (options.isValid()) {

			AssemblerSettings assemblerSettings = new AssemblerSettings();

			assemblerSettings.setKmerSize(options.getKmerSize());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings
					.setMinEdgeFrequency(options.getMinEdgeFrequency());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMinEdgeRatio(options.getMinEdgeRatio());
			assemblerSettings.setMaxPotentialContigs(options
					.getMaxPotentialContigs());
			assemblerSettings.setMinContigRatio(options.getMinContigRatio());
			assemblerSettings.setMinUniqueReads(options.getMinUniqueReads());

			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setReferenceDir(options.getReferenceDir());
			realigner.setRegionsGtf(options.getTargetRegionFile());
			realigner.setTempDir(options.getWorkingDir());
			realigner.setAssemblerSettings(assemblerSettings);
			realigner.setNumThreads(options.getNumThreads());
			realigner.setMinContigMapq(options.getMinContigMapq());
			realigner.setAllowedMismatchesFromContig(options.getAllowedMismatchesFromContig());
			realigner.setShouldReprocessUnaligned(!options.isSkipUnalignedAssembly());

			long s = System.currentTimeMillis();

			realigner.reAlign(options.getInputFile(), options.getOutputFile());

			long e = System.currentTimeMillis();

			System.out.println("Elapsed seconds: " + (e - s) / 1000);
		}
	}

	public static void main(String[] args) throws Exception {
		ReAligner realigner = new ReAligner();

		long s = System.currentTimeMillis();

/*
		String input = "/home/lisle/ayc/sim/sim261/chr1/sorted.bam";
		String output = "/home/lisle/ayc/sim/sim261/chr1/realigned.bam";
		String reference = "/home/lisle/reference/chr1/chr1.fa";
		String regions = "/home/lisle/ayc/regions/chr1_261.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim261/chr1/working";
*/
		
/*
		String input = "/home/lisle/ayc/sim/sim261/chr17/sorted.bam";
		String output = "/home/lisle/ayc/sim/sim261/chr17/realigned.bam";
		String reference = "/home/lisle/reference/chr17/chr17.fa";
		String regions = "/home/lisle/ayc/regions/chr17_261.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim261/chr17/working";
*/		
		String input = "/home/lmose/dev/ayc/sim/sim261/chr11/sorted.bam";
		String output = "/home/lmose/dev/ayc/sim/sim261/chr11/realigned.bam";
		String reference = "/home/lmose/reference/chr11/chr11.fa";
		String regions = "/home/lmose/dev/ayc/regions/chr11_261.gtf";
		String tempDir = "/home/lmose/dev/ayc/sim/sim261/chr11/working";


		
/*
		String input = "/home/lisle/ayc/sim/sim261/chr13/sorted.bam";
		String output = "/home/lisle/ayc/sim/sim261/chr13/realigned.bam";
		String reference = "/home/lisle/reference/chr13/chr13.fa";
		String regions = "/home/lisle/ayc/regions/chr13_261.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim261/chr13/working";
*/
		
/*		
		String input = "/home/lisle/ayc/sim/sim261/chr16/sorted.bam";
		String output = "/home/lisle/ayc/sim/sim261/chr16/realigned.bam";
		String reference = "/home/lisle/reference/chr16/chr16.fa";
		String regions = "/home/lisle/ayc/regions/chr16_261.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim261/chr16/working";
*/
		
/*
		String input = "/home/lisle/ayc/sim/sim261/sorted.bam";
		String output = "/home/lisle/ayc/sim/sim261/realigned.bam";
		String reference = "/home/lisle/reference/chr13/chr13.fa";
		String regions = "/home/lisle/ayc/regions/chr13_261.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim261/working";
*/
/*		
		String input = "/home/lisle/ayc/sim/sim1/bug/chr8_141889351_141889791.bam";
		String output = "/home/lisle/ayc/sim/sim1/bug/realigned.bam";
		String reference = "/home/lisle/reference/chr8/chr8.fa";
		String regions = "/home/lisle/ayc/regions/chr8_141889351_141889791.gtf";
		String tempDir = "/home/lisle/ayc/sim/sim1/bug/working";
*/


		AssemblerSettings settings = new AssemblerSettings();
		settings.setKmerSize(63);
		settings.setMinContigLength(100);
		settings.setMinEdgeFrequency(3);
		settings.setMinNodeFrequncy(3);
		settings.setMinEdgeRatio(.02);
		settings.setMaxPotentialContigs(10000);
		settings.setMinContigRatio(.3);
		settings.setMinUniqueReads(1);

		realigner.setAssemblerSettings(settings);
		
		realigner.setMinContigMapq(1);
		realigner.setAllowedMismatchesFromContig(2);
		realigner.setReference(reference);
		realigner.setRegionsGtf(regions);
		realigner.setTempDir(tempDir);
		realigner.setNumThreads(1);

		realigner.reAlign(input, output);

		long e = System.currentTimeMillis();

		System.out.println("Elapsed seconds: " + (e - s) / 1000);
	}
}
