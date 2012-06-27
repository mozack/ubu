package edu.unc.bioinf.ubu.assembly;

public class Contig {

	private String descriptor;
	private String sequence;
	
	public Contig(String descriptor, String sequence) {
		this.descriptor = descriptor;
		this.sequence = sequence;
	}
	
	public String getDescriptor() {
		return descriptor;
	}
	
	public String getSequence() {
		return sequence;
	}
}
