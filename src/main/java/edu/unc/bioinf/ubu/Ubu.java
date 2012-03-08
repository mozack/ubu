package edu.unc.bioinf.ubu;

import java.util.Arrays;

import edu.unc.bioinf.ubu.fastq.FastqMapsplicePrep;
import edu.unc.bioinf.ubu.sam.GenomeToTranscriptome;

/**
 * Entry point for Ubu java command line utilities
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Ubu {
	
	private static final String SAM_DIFF = "samdiff";
	private static final String SAM_FILTER = "samfilter";
	private static final String JUNC = "junc";
	private static final String MAPSPLICE_PREP = "mapsplice-prep";
	
	public void run(String[] args) throws Exception {
		if (args.length == 0) {
			printAvailablePrograms();
		} else {
		
			String cmd = args[0];
			
			String[] argz = Arrays.copyOfRange(args, 1, args.length);
			
			if (cmd.equals("xlate")) {
				GenomeToTranscriptome.run(argz);
			} else if (cmd.equals(SAM_DIFF)) {
				printNotYetSupported(cmd);
			} else if (cmd.equals(SAM_FILTER)) {
				printNotYetSupported(cmd);
			} else if (cmd.equals(JUNC)) {
				printNotYetSupported(cmd);
			} else if (cmd.equals(MAPSPLICE_PREP)) {
				FastqMapsplicePrep.run(argz);
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
		
		printProgram("xlate", "\tTranslate from genome to transcriptome coordinates");
		printProgram("samdiff", "\tDiff two SAM/BAM files outputting discrepancies in corresponding SAM/BAM files");
		printProgram("samfilter", "Filter reads from a SAM or BAM file");
		printProgram("junc", "\tCount splice junctions in a SAM or BAM file");
		printProgram("mapsplice-prep", "Prep a single fastq file for Mapsplice (for Casava 1.8 output)");
	}
	
	private void printProgram(String name, String desc) {
		System.out.println("\t" + name + "\t" + desc);
	}
	
	public static void main(String[] args) throws Exception {
		new Ubu().run(args);
	}
}
