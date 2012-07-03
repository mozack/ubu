package edu.unc.bioinf.ubu.assembly;

import net.sf.samtools.SAMRecord;

public class ReadPosition {

	private SAMRecord read;
	private int position;
	
	public ReadPosition(SAMRecord read, int position) {
		this.read = read;
		this.position = position;
	}
	
	public SAMRecord getRead() {
		return read;
	}
	
	public int getPosition() {
		return position;
	}
}
