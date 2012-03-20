package edu.unc.bioinf.ubu;

import java.util.Arrays;

import edu.unc.bioinf.ubu.fastq.FastqFormatter;
import edu.unc.bioinf.ubu.sam.GenomeToTranscriptome;
import edu.unc.bioinf.ubu.sam.SamFileDiff;
import edu.unc.bioinf.ubu.sam.SamSummarizer;

/**
 * Entry point for Ubu java command line utilities
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Ubu {
	
	private static final String TRANSLATE = "xlate";
	private static final String SAM_DIFF = "sam-diff";
	private static final String SAM_FILTER = "sam-filter";
	private static final String SAM_SUMMARIZE = "sam-summary";
	private static final String JUNC = "junc";
	private static final String FASTQ_FORMAT = "fastq-format";
	
	private static final int MAX_CMD_LEN = 15;
	
	public void run(String[] args) throws Exception {
		if (args.length == 0) {
			printAvailablePrograms();
		} else {
		
			String cmd = args[0];
			
			String[] argz = Arrays.copyOfRange(args, 1, args.length);
			
			if (cmd.equals(TRANSLATE)) {
				GenomeToTranscriptome.run(argz);
			} else if (cmd.equals(SAM_DIFF)) {
				SamFileDiff.run(argz);
			} else if (cmd.equals(SAM_FILTER)) {
				printNotYetSupported(cmd);
			} else if (cmd.equals(JUNC)) {
				printNotYetSupported(cmd);
			} else if (cmd.equals(FASTQ_FORMAT)) {
				FastqFormatter.run(argz);
			} else if (cmd.equals(SAM_SUMMARIZE)) {
				SamSummarizer.run(argz);
			} else {
				System.out.println("Command [" + cmd + "] is unrecognized.");
				printAvailablePrograms();
			}
		}
	}
	
	private void printNotYetSupported(String cmd) {
		System.out.println("Sorry, [" + cmd + "] is not yet supported.  Check back soon.");
	}
	
	private void printAvailablePrograms() {
		System.out.println("UNC-Chapel Hill Bioinformatics Utilities");
		System.out.println("Available commands:");
		
		printProgram(getPaddedString(TRANSLATE), "Translate from genome to transcriptome coordinates");
		printProgram(getPaddedString(SAM_DIFF), "Diff two SAM/BAM files outputting discrepancies in corresponding SAM/BAM files");
		printProgram(getPaddedString(SAM_FILTER), "Filter reads from a SAM or BAM file");
		printProgram(getPaddedString(SAM_SUMMARIZE), "Output summary statistics per reference for a SAM/BAM file.");
		printProgram(getPaddedString(JUNC), "Count splice junctions in a SAM or BAM file");
		printProgram(getPaddedString(FASTQ_FORMAT), "Format a single FASTQ file (clean up read ids and/or convert quality scoring)");
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
