package edu.unc.bioinf.ubu.assembly;

public class AssemblerSettings {

	private int kmerSize;
	private int minEdgeFrequency;
	private int minNodeFrequncy;
	private int minContigLength;
	private double minEdgeRatio;
	
	public int getKmerSize() {
		return kmerSize;
	}
	public void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
	}
	public int getMinEdgeFrequency() {
		return minEdgeFrequency;
	}
	public void setMinEdgeFrequency(int minEdgeFrequency) {
		this.minEdgeFrequency = minEdgeFrequency;
	}
	public int getMinNodeFrequncy() {
		return minNodeFrequncy;
	}
	public void setMinNodeFrequncy(int minNodeFrequncy) {
		this.minNodeFrequncy = minNodeFrequncy;
	}
	public int getMinContigLength() {
		return minContigLength;
	}
	public void setMinContigLength(int minContigLength) {
		this.minContigLength = minContigLength;
	}
	public double getMinEdgeRatio() {
		return minEdgeRatio;
	}
	public void setMinEdgeRatio(double minEdgeRatio) {
		this.minEdgeRatio = minEdgeRatio;
	}
	
	public String getDescription() {
		StringBuffer str = new StringBuffer();
		
		appendSetting(str, "kmerSize", kmerSize);
		appendSetting(str, "minEdgeFrequency", minEdgeFrequency);
		appendSetting(str, "minNodeFrequncy", minNodeFrequncy);
		appendSetting(str, "minContigLength", minContigLength);
		appendSetting(str, "minEdgeRatio", minEdgeRatio);
		
		return str.toString();
	}
	
	private void appendSetting(StringBuffer str, String setting, int value) {
		str.append(setting);
		str.append(": ");
		str.append(value);
		str.append('\n');
	}
	
	private void appendSetting(StringBuffer str, String setting, double value) {
		str.append(setting);
		str.append(": ");
		str.append(value);
		str.append('\n');
	}
}
