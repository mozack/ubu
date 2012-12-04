package edu.unc.bioinf.ubu;

import java.util.Arrays;

import edu.unc.bioinf.ubu.assembly.ReAligner;
import edu.unc.bioinf.ubu.fastq.FastqFormatter;
import edu.unc.bioinf.ubu.fastq.Sam2Fastq;
import edu.unc.bioinf.ubu.sam.GenomeToTranscriptome;
import edu.unc.bioinf.ubu.sam.SAMFilter;
import edu.unc.bioinf.ubu.sam.SamConverter;
import edu.unc.bioinf.ubu.sam.SamFileDiff;
import edu.unc.bioinf.ubu.sam.SamSummarizer;
import edu.unc.bioinf.ubu.sam.SpliceJunctionCounter;

/**
 * Entry point for Ubu java command line utilities
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Ubu {
	
	private static final String TRANSLATE = "sam-xlate";
	private static final String SAM_DIFF = "sam-diff";
	private static final String SAM_FILTER = "sam-filter";
	private static final String SAM_SUMMARIZE = "sam-summary";
	private static final String SAM_CONVERT = "sam-convert";
	private static final String JUNC = "sam-junc";
	private static final String SAM2FASTQ = "sam2fastq";
	private static final String FASTQ_FORMAT = "fastq-format";
	private static final String REALIGN = "realign";
	
	private static final int MAX_CMD_LEN = 15;
	
	public void run(String[] args) throws Exception {
		if (args.length == 0) {
			printAvailablePrograms();
		} else {
			
			System.out.println("Java version: " + System.getProperty("java.version"));
			System.out.println("Max memory: " + Runtime.getRuntime().maxMemory());
			System.out.println("Total memory: " + Runtime.getRuntime().totalMemory());
			System.out.println("Free memory: " + Runtime.getRuntime().freeMemory());
		
			String cmd = args[0];
			
			String[] argz = Arrays.copyOfRange(args, 1, args.length);
			
			if (cmd.equals(TRANSLATE)) {
				GenomeToTranscriptome.run(argz);
			} else if (cmd.equals(SAM_DIFF)) {
				SamFileDiff.run(argz);
			} else if (cmd.equals(SAM_FILTER)) {
				SAMFilter.run(argz);
			} else if (cmd.equals(JUNC)) {
				SpliceJunctionCounter.run(argz);
			} else if (cmd.equals(FASTQ_FORMAT)) {
				FastqFormatter.run(argz);
			} else if (cmd.equals(SAM_SUMMARIZE)) {
				SamSummarizer.run(argz);
			} else if (cmd.equals(SAM_CONVERT)) {
				SamConverter.run(argz);
			} else if (cmd.equals(SAM2FASTQ)) {
				Sam2Fastq.run(argz);
			} else if (cmd.equals(REALIGN)) {
				ReAligner.run(argz);
			} else {
				System.out.println("Command [" + cmd + "] is unrecognized.");
				printAvailablePrograms();
			}
		}
	}
	
	private void printAvailablePrograms() {
		System.out.println("UNC-Chapel Hill Bioinformatics Utilities");
		System.out.println("Available commands:");
		
		printProgram(getPaddedString(TRANSLATE), "Translate from genome to transcriptome coordinates");
		printProgram(getPaddedString(SAM_DIFF), "Diff two SAM/BAM files outputting discrepant reads in corresponding SAM/BAM files");
		printProgram(getPaddedString(SAM_FILTER), "Filter reads from a paired end SAM or BAM file (only outputs paired reads)");
		printProgram(getPaddedString(SAM_SUMMARIZE), "Output summary statistics per reference for a SAM/BAM file (Aligned reads only).");
		printProgram(getPaddedString(SAM_CONVERT), "Convert SAM/BAM file content (i.e. convert quality from phred64 to phred33)");
		printProgram(getPaddedString(JUNC), "Count splice junctions in a SAM or BAM file");
		printProgram(getPaddedString(SAM2FASTQ), "Convert SAM/BAM file to FASTQ");
		printProgram(getPaddedString(FASTQ_FORMAT), "Format a single FASTQ file (clean up read ids and/or convert quality scoring)");
		printProgram(getPaddedString(REALIGN), "");
	}
	
	private String getPaddedString(String str) {
		int pad = MAX_CMD_LEN - str.length();
		for (int i=0; i<pad; i++) {
			str += ' ';
		}
		return str;
	}
	
	private void printProgram(String name, String desc) {
		System.out.println("\t" + name + "\t" + desc);
	}
	
	public static void main(String[] args) throws Exception {
		new Ubu().run(args);
	}
}
