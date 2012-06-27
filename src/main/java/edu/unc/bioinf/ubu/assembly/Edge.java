package edu.unc.bioinf.ubu.assembly;

public class Edge {

	private Node from;
	private Node to;
	private int count = 1;
	private int hashCode;
	
	public Edge(Node from, Node to) {
		this.from = from;
		this.to = to;
		hashCode = from.hashCode() * 37 + to.hashCode();
	}
	
	public String toString() {
		return count + "_" + from.toString() + "->" + to.toString();
	}
	
	public boolean equals(Object object) {
		Edge that = (Edge) object;
		return this.from.equals(that.from) && this.to.equals(that.to);
	}
	
	public int hashCode() {
		return hashCode;
	}
	
	public void incrementCount() {
		count++;
	}
	
	public int getCount() {
		return count;
	}
	
	public Node getTo() {
		return to;
	}
	
	public Node getFrom() {
		return from;
	}
	
	public void remove() {
		to.removeFromEdge(this);
		from.removeToEdge(this);
	}
}
