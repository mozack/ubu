package edu.unc.bioinf.ubu.assembly;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Node {

	private Sequence sequence;
	private int count = 1;
	private Map<Node, Edge> toEdges = new HashMap<Node, Edge>(2);
	private Map<Node, Edge> fromEdges = new HashMap<Node, Edge>(2);
	private String contributingRead = null;
	private boolean hasMultipleUniqueReads = false;
	
	public Node(String sequence) {
		this(new Sequence(sequence));
	}
	
	public Node(Sequence sequence) {
		this.sequence = sequence;
	}
	
	public String toString() {
		return count + "_" + sequence;
	}
	
	public int hashCode() {
		return sequence.hashCode();
	}
	
	public boolean equals(Object object) {
		Node that = (Node) object;
		return this.sequence.equals(that.sequence);
	}
	
	public void incrementCount() {
		count++;
	}
	
	public int getCount() {
		return count;
	}
	
	public Collection<Edge> getToEdges() {
		return toEdges.values();
	}
	
	public Collection<Edge> getFromEdges() {
		return fromEdges.values();
	}
	
	public Sequence getSequence() {
		return sequence;
	}
	
	public void addReadSequence(String sequence) {
		if (!hasMultipleUniqueReads) {
			if (contributingRead == null) {
				contributingRead = sequence;
			} else {
				if (!contributingRead.equals(sequence)) {
					hasMultipleUniqueReads = true;
					contributingRead = null;
				}
			}
		}
	}
	
	private void printMultiEdges() {
		int aboveThreshold = 0;
		
		for (Edge edge : toEdges.values()) {
			if (edge.getCount() > 50) {
				aboveThreshold++;
			}
		}
		
		if (aboveThreshold > 1) {
			System.out.println("------- Edge --------");
			for (Edge edge : toEdges.values()) {
				System.out.print(edge.getCount() + ", ");
			}
			
			System.out.println();
		}
	}
	
	public List<Edge> getInfrequentEdges(double minFreq) {
		List<Edge> infrequentEdges = new ArrayList<Edge>();
		
		infrequentEdges.addAll(getInfrequentEdges(minFreq, toEdges.values()));
		infrequentEdges.addAll(getInfrequentEdges(minFreq, fromEdges.values()));
		
		return infrequentEdges;
	}
	
	private List<Edge> getInfrequentEdges(double minFreq, Collection<Edge> edges) {
		List<Edge> infrequentEdges = new ArrayList<Edge>();
		
		double total = getEdgeTotal(edges);
				
		for (Edge edge : edges) {
			if (((double) edge.getCount() / total) < minFreq) {
				infrequentEdges.add(edge);
			}
		}
		
		return infrequentEdges;
	}
	
	private double getEdgeTotal(Collection<Edge> edges) {
		double total = 0.0;
		
		for (Edge edge : edges) {
			total += edge.getCount();
		}

		return total;
	}
	
	public List<Edge> getFrequentToEdges(double minFreq) {
		
		List<Edge> frequentEdges = new ArrayList<Edge>();
		
		double total = getEdgeTotal(toEdges.values());
		
		for (Edge edge : toEdges.values()) {
			if (((double) edge.getCount() / total) >= minFreq) {
				frequentEdges.add(edge);
			}
		}
		
		return frequentEdges;
	}
	
	public boolean hasMultipleUniqueReads() {
		return hasMultipleUniqueReads;
	}
	
	public Edge getMostCommonEdge() {
		printMultiEdges();
		
		Edge topEdge = null;
		int freq = 0;
		
//		if (toEdges.size() > 1)
//			System.out.println("------- Edge --------");
		
		for (Edge edge : toEdges.values()) {
//			if (toEdges.size() > 1) 
//				System.out.print(edge.getCount() + ", ");
			
			if (edge.getCount() > freq) {
				topEdge = edge;
				freq = edge.getCount();
			}
		}
		
//		if (toEdges.size() > 1)
//			System.out.println();
		
		return topEdge;
	}
	
	public void addToEdge(Node to) {
		Edge edge = toEdges.get(to);
		
		if (edge == null) {
			edge = new Edge(this, to);
			toEdges.put(to, edge);
			to.updateFromEdges(edge);
		} else {
			edge.incrementCount();
		}
	}
	
	public boolean isRootNode() {
		return fromEdges.isEmpty();
	}
	
	private void updateFromEdges(Edge edge) {
		fromEdges.put(edge.getFrom(), edge);
	}
	
	public boolean isSingleton() {
		return fromEdges.isEmpty() && toEdges.isEmpty();
	}
	
	public void removeToEdge(Edge edge) {
		this.toEdges.remove(edge.getTo());
	}
	
	public void removeFromEdge(Edge edge) {
		this.fromEdges.remove(edge.getFrom());
	}
}
