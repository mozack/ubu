package edu.unc.bioinf.ubu.assembly;

import java.io.IOException;

public class Aligner {
	
	private String reference;
	
	public Aligner(String reference) {
		this.reference = reference;
	}

	public void align(String input, String outputSam) throws IOException, InterruptedException {
		String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.out.println("Running: [" + cmd + "]");
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		if (ret != 0) {
			throw new RuntimeException("BWA exited with non-zero return code : " + ret);
		}
	}
}
