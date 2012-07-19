package edu.unc.bioinf.ubu.assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringWriter;
import java.util.ArrayList;
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
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import edu.unc.bioinf.ubu.gtf.Feature;
import edu.unc.bioinf.ubu.gtf.GtfLoader;
import edu.unc.bioinf.ubu.sam.ReadBlock;

public class ReAligner {
	
	private Assembler assembler;
	
	private Aligner aligner;
	
	//private Aligner aligner = new Aligner("/home/lisle/reference/chr5/chr5.fa");
	//private Aligner aligner = new Aligner("/home/lisle/reference/chrX/chrX.fa");
	
	private Set<SAMRecord> updatedReads = new HashSet<SAMRecord>();
	
	private SAMFileHeader samHeader;
	
	private SAMFileWriter outputReadsBam;
	
	private long startMillis;
	
	private List<Feature> regions;
	
	private String regionsGtf;
	
	private String tempDir;
	
	private String reference;
	
	private AssemblerSettings assemblerSettings;
	
	public void reAlign(String inputSam, String outputSam) throws Exception {
		
		System.out.println("input: " + inputSam);
		System.out.println("output: " + outputSam);
		System.out.println("regions: " + regionsGtf);
		System.out.println("reference: " + reference);
		System.out.println("working dir: " + tempDir);
		System.out.println(assemblerSettings.getDescription());
		
		startMillis = System.currentTimeMillis();
		
		//String contigsFasta = outputPrefix + "_contigs.fasta";
		//String contigsSam = outputPrefix + "_contigs.sam";
		
		init();
		
		log("Loading target regions");
		loadRegions();
		
		log("Reading Input SAM Header");
		getSamHeader(inputSam);
		
		log("Initializing output SAM File");
		initOutputFile(outputSam);
		
		for (Feature region : regions) {	
			
			log("Extracting targeted region: " + region.getDescriptor());
			String targetRegionBam = extractTargetRegion(inputSam, region);
			
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			String contigsSam   = tempDir + "/" + region.getDescriptor() + "_contigs.sam";
			
			log("Initializing assembler");
			initAssembler();
			
			log("Assembling contigs");
			List<Contig> contigs = assembler.assembleContigs(targetRegionBam, contigsFasta);
			
			if (contigs.size() > 0) {
				log("Aligning contigs");
				aligner.align(contigsFasta, contigsSam);
				
				log("Adjusting reads");
				adjustReads(contigsSam, contigs);
				
				log("Writing adjusted reads");
				outputReads(outputSam);
			} else {
				log ("No contigs assembled for region: " + region.getDescriptor());
			}
			
			updatedReads.clear();
		}
		
		log("Closing output BAM");
		outputReadsBam.close();
		
		System.out.println("Done.");
	}
	
	private String extractTargetRegion(String inputSam, Feature region) throws IOException, InterruptedException {
		String extractFile = tempDir + "/" + region.getDescriptor() + ".bam";
		
		String location = region.getSeqname() + ":" + region.getStart() + "-" + region.getEnd();
		
		String cmd = "samtools view -b " + inputSam + " " + location + " -o " + extractFile;
		
		System.out.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);

		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.out.println("Extract time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			String stdout = getOutput(proc.getInputStream());
			String stderr = getOutput(proc.getErrorStream());
			System.out.println("Samtools stdout: " + stdout);
			System.out.println("Samtools stderr: " + stderr);
			
			throw new RuntimeException("Samtools exited with non-zero return code : " + ret);
		}
		
		return extractFile;
	}
	
	private String getOutput(InputStream is) throws IOException {
		StringWriter writer = new StringWriter();
		
		Reader reader = new BufferedReader(
                new InputStreamReader(is));
		
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

        reader.close();
	}
	
	private void adjustReads(String contigSam, List<Contig> contigs) {
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
        		//TODO: Handle multi-mappers (update XH tags?)
    			SAMRecord updatedRead = updateReadAlignment(contigReadBlocks, readPosition);
    			if (updatedRead != null) {
    				updatedReads.add(updatedRead);
    			}
        	}
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
		
		if (blocks.size() > 0) {
			int newAlignmentStart = blocks.get(0).getReferenceStart();
			String newCigar = ReadBlock.toCigarString(blocks);
			
			read.setCigarString(newCigar);
			read.setAlignmentStart(newAlignmentStart);
		} else {
			//TODO: Investigate how this could happen.
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
		
        outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samHeader,
                true, new File(outputReadsBamFilename));
    }
    
	private void outputReads(String readsBam) {		
		System.out.println("Writing " + updatedReads.size() + " reads.");
        for (SAMRecord read : updatedReads) {
        	outputReadsBam.addAlignment(read);
        }
	}
	
	private void initAssembler() {
		assembler = new Assembler();
		
		assembler.setKmerSize(assemblerSettings.getKmerSize());
		assembler.setMinEdgeFrequency(assemblerSettings.getMinEdgeFrequency());
		assembler.setMinNodeFrequncy(assemblerSettings.getMinNodeFrequncy());
		assembler.setMinContigLength(assemblerSettings.getMinContigLength());
		assembler.setMinEdgeRatio(assemblerSettings.getMinEdgeRatio());
	}
	
	private void init() {
//		reference = "/home/lisle/reference/chr17/chr17.fa";
//		regionsGtf = "/home/lisle/ayc/regions/chr17.gtf";
//		tempDir = "/home/lisle/ayc/case0/round2/working1";
		
		aligner = new Aligner(reference);
		
		File workingDir = new File(tempDir);
		if (workingDir.exists()) {
			if (!workingDir.delete()) {
				throw new IllegalStateException("Unable to delete: " + tempDir);
			}
		}
		
		if (!workingDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir);
		}
	}
	
	public void setReference(String reference) {
		this.reference = reference;
	}
	
	public void setTempDir(String temp) {
		this.tempDir = temp;
	}
	
	public void setAssemblerSettings(AssemblerSettings settings) {
		this.assemblerSettings = settings;
	}
	
	public static void run(String[] args) throws Exception {
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			
			AssemblerSettings assemblerSettings = new AssemblerSettings();
			
			assemblerSettings.setKmerSize(options.getKmerSize());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings.setMinEdgeFrequency(options.getMinEdgeFrequency());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMinEdgeRatio(options.getMinEdgeRatio());
			
			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setRegionsGtf(options.getTargetRegionFile());
			realigner.setTempDir(options.getWorkingDir());
			realigner.setAssemblerSettings(assemblerSettings);
			
			long s = System.currentTimeMillis();
			
			realigner.reAlign(options.getInputFile(), options.getOutputFile());
			
			long e = System.currentTimeMillis();
			
			System.out.println("Elapsed seconds: " + (e-s)/1000);
		}
	}
	
	public static void main(String[] args) throws Exception {
		ReAligner realigner = new ReAligner();
		
		long s = System.currentTimeMillis();
		
		String input     = "/home/lisle/ayc/case4/tumor/chr15_99503528_99507968.bam";
		String output    = "/home/lisle/ayc/case4/tumor/realigned.bam";
		String reference = "/home/lisle/reference/chr15/chr15.fa";
		String regions   = "/home/lisle/ayc/regions/chr15_99503528_99507968.gtf";
		String tempDir   = "/home/lisle/ayc/case4/tumor/working";
		
		AssemblerSettings settings = new AssemblerSettings();
		settings.setKmerSize(47);
		settings.setMinContigLength(101);
		settings.setMinEdgeFrequency(10);
		settings.setMinNodeFrequncy(10);
		settings.setMinEdgeRatio(.02);
		
		realigner.setAssemblerSettings(settings);
		
//		reference = "/home/lisle/reference/chr17/chr17.fa";
//		regionsGtf = "/home/lisle/ayc/regions/chr17.gtf";
//		tempDir = "/home/lisle/ayc/case0/round2/working1";

		
//		realigner.reAlign("/home/lisle/ayc/case0/round2/case0_tumor.bam", "/home/lisle/ayc/case0/round2/full.bam");
		
		realigner.setReference(reference);
		realigner.setRegionsGtf(regions);
		realigner.setTempDir(tempDir);
		
		realigner.reAlign(input, output);
		
		long e = System.currentTimeMillis();
		
		System.out.println("Elapsed seconds: " + (e-s)/1000);
	}
}
