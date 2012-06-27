package edu.unc.bioinf.ubu.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Counts implements Cloneable {

	private int totalEdgeCount = 0;
	private List<Integer> edgeCounts = new ArrayList<Integer>();
	private boolean isTerminatedAtRepeat = false;
	private boolean isSorted = false;
	
	public void incrementEdgeCounts(int count) {
		totalEdgeCount += count;
		edgeCounts.add(count);
		isSorted = false;
	}
	
	public int getTotalEdgeCount() {
		return totalEdgeCount;
	}
	
	public int getMedianEdgeCount() {
		if (edgeCounts.size() == 0) {
			return 0;
		}
		
		if (!isSorted) {
			Collections.sort(edgeCounts);
		}
		
		return edgeCounts.get(edgeCounts.size()/2);
	}
	
	public int getMinEdgeCount() {
		if (edgeCounts.size() == 0) {
			return 0;
		}
		if (!isSorted) {
			Collections.sort(edgeCounts);
		}
		
		return edgeCounts.get(0);
	}
	
	public int getNumEdges() {
		return edgeCounts.size();
	}
	
	public void setTerminatedAtRepeat(boolean isTerminatedAtRepeat) {
		this.isTerminatedAtRepeat = isTerminatedAtRepeat;
	}
	
	public boolean isTerminatedAtRepeat() {
		return isTerminatedAtRepeat;
	}
	
	public Object clone() {
		Counts clone = new Counts();
		clone.edgeCounts.addAll(this.edgeCounts);
		clone.totalEdgeCount = this.totalEdgeCount;
		
		return clone;
	}
	
	public String toString() {
		return 
			"_numedges:" + getNumEdges() + 
			"_totaledgecounts:" + getTotalEdgeCount() +
			"_medianedgecount:" + getMedianEdgeCount() +
			"_minedgecount:" + getMinEdgeCount() + 
			"_terminatedatrepeat:" + isTerminatedAtRepeat; 
	}
}