package edu.unc.bioinf.ubu.assembly;

import java.io.IOException;

public class Aligner {
	
	private String reference;
	private int numThreads;
	
	public Aligner(String reference, int numThreads) {
		this.reference = reference;
		this.numThreads = numThreads;
	}
	
	public void align(String input, String outputSam) throws IOException, InterruptedException {
		String cmd = "bwa bwasw -t " + numThreads + " -f " + outputSam + " " + reference + " " + input;
		
		runCommand(cmd);
	}

	private void runCommand(String cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.out.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.out.println("BWA time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("BWA exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
	}
	
	public void shortAlign(String input, String outputSam) throws IOException, InterruptedException {
		String sai = outputSam + ".sai";
		
		String aln = "bwa aln " + reference + " " + input + " -f " + sai + " -t " + numThreads;
		
		runCommand(aln);
		
		String convert = "bwa samse " + reference + " " + sai + " " + input + " -f " + outputSam  + " -n 100";
		
		runCommand(convert);
	}
	
	public void index() throws IOException, InterruptedException {
		runCommand("bwa index " + reference);
	}
}
