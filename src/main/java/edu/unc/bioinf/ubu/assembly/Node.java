package edu.unc.bioinf.ubu.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Node {

	private Sequence sequence;
	private int count = 1;
	private Edge[] toEdges = new Edge[0];
	private Edge[] fromEdges = new Edge[0];
//	private List<Edge> toEdges = new ArrayList<Edge>(1);
//	private List<Edge> fromEdges = new ArrayList<Edge>(1);
//	private Map<Node, Edge> toEdges = new HashMap<Node, Edge>(2);
//	private Map<Node, Edge> fromEdges = new HashMap<Node, Edge>(2);
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
		return Arrays.asList(toEdges);
	}
	
	public Collection<Edge> getFromEdges() {
		return Arrays.asList(fromEdges);
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
		
		for (Edge edge : toEdges) {
			if (edge.getCount() > 50) {
				aboveThreshold++;
			}
		}
		
		if (aboveThreshold > 1) {
			System.out.println("------- Edge --------");
			for (Edge edge : toEdges) {
				System.out.print(edge.getCount() + ", ");
			}
			
			System.out.println();
		}
	}
	
	public List<Edge> getInfrequentEdges(double minFreq) {
		List<Edge> infrequentEdges = new ArrayList<Edge>();
		
		infrequentEdges.addAll(getInfrequentEdges(minFreq, toEdges));
		infrequentEdges.addAll(getInfrequentEdges(minFreq, fromEdges));
		
		return infrequentEdges;
	}
	
	private List<Edge> getInfrequentEdges(double minFreq, Edge[] edges) {
		List<Edge> infrequentEdges = new ArrayList<Edge>();
		
		double total = getEdgeTotal(edges);
				
		for (Edge edge : edges) {
			if (((double) edge.getCount() / total) < minFreq) {
				infrequentEdges.add(edge);
			}
		}
		
		return infrequentEdges;
	}
	
	private double getEdgeTotal(Edge[] edges) {
		double total = 0.0;
		
		for (Edge edge : edges) {
			total += edge.getCount();
		}

		return total;
	}
	
	public List<Edge> getFrequentToEdges(double minFreq) {
		
		List<Edge> frequentEdges = new ArrayList<Edge>();
		
		double total = getEdgeTotal(toEdges);
		
		for (Edge edge : toEdges) {
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
		
		for (Edge edge : toEdges) {
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
		
	private Edge[] addEdge(Edge[] edges, Edge edge) {
//		int minCapacity = edges.length + 1;
		//Edge[] newEdges = new Edge[edges.length + 1];
		Edge[] newEdges = Arrays.copyOf(edges, edges.length + 1);
		newEdges[edges.length] = edge;
		
		return newEdges;
		
//        int oldCapacity = edges.length;
//        if (minCapacity > oldCapacity) {
//            Object oldData[] = edges;
//            int newCapacity = (oldCapacity * 3)/2 + 1;
//            if (newCapacity < minCapacity)
//                newCapacity = minCapacity;
//            // minCapacity is usually close to size, so this is a win:
//            elementData = Arrays.copyOf(elementData, newCapacity);
//        }
	}
	
	public void addToEdge(Node to) {
		Edge edge = findEdge(to);
		
		if (edge == null) {
			edge = new Edge(this, to);
			toEdges = addEdge(toEdges, edge);
			to.updateFromEdges(edge);
		} else {
			edge.incrementCount();
		}
	}
	
	public boolean isRootNode() {
		return fromEdges.length == 0;
	}
	
	private void updateFromEdges(Edge edge) {
		fromEdges = addEdge(fromEdges, edge);
	}
	
	public boolean isSingleton() {
		return fromEdges.length == 0 && toEdges.length == 0;
	}
	
	//TODO: Smarter array allocation
	public void removeToEdge(Edge edge) {
//		this.toEdges.remove(edge.getTo());
		int index = this.findToEdgeIndex(edge.getTo());
		if (index > -1) {
			Edge[] newEdges = new Edge[toEdges.length-1];
			System.arraycopy(toEdges, 0, newEdges, 0, index);
			if ((index) < newEdges.length) { 
				System.arraycopy(toEdges, index+1, newEdges, index, newEdges.length-index);
			}
			
			toEdges = newEdges;
		}
	}
	
	public void removeFromEdge(Edge edge) {
//		this.fromEdges.remove(edge.getFrom());
		int index = this.findFromEdgeIndex(edge.getFrom());
		if (index > -1) {
			Edge[] newEdges = new Edge[fromEdges.length-1];
			System.arraycopy(fromEdges, 0, newEdges, 0, index);
			if ((index) < newEdges.length) {
				System.arraycopy(fromEdges, index+1, newEdges, index, newEdges.length-index);
			}
			
			fromEdges = newEdges;
		}
	}
	
	private int findToEdgeIndex(Node to) {
		for (int i=0; i<toEdges.length; i++) {
			if (toEdges[i].getTo().equals(to)) {
				return i;
			}
		}
		
		return -1;
	}
	
	private int findFromEdgeIndex(Node from) {
		for (int i=0; i<fromEdges.length; i++) {
			if (fromEdges[i].getFrom().equals(from)) {
				return i;
			}
		}
		
		return -1;
	}

	
	private Edge findEdge(Node to) {
		for (Edge edge : toEdges) {
			if (edge.getTo().equals(to)) {
				return edge;
			}
		}
		
		return null;
	}
}
