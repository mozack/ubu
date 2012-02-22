package edu.unc.bioinf.ubu.sam;

/**
 * Represents a set of coordinates specified by a start and stop position.
 * The coordinates are assumed to be inclusive.
 * 
 * @author lmose
 */
public class Coordinate implements Comparable<Coordinate> {

    private final int start;
    private final int stop;
    
    public Coordinate(int start, int stop) {
        this.start = start;
        this.stop = stop;
    }
    
    public int getStart() {
        return start;
    }
    
    public int getStop() {
        return stop;
    }
    
    public boolean contains(int position) {
        return (position >= start) && (position <= stop);
    }
    
    public int getLength() {
        return stop-start+1;
    }
    
    public String toString() {
        return start + "," + stop;
    }
    
    @Override
    public boolean equals(Object object) {
        Coordinate that = (Coordinate) object;
        return this.start == that.start && this.stop == that.stop;
    }
    
    @Override
    public int hashCode() {
        int result = 17;
        
        result = 31 * result + start;
        result = 31 * result + stop;
        
        return result;
    }

    @Override
    public int compareTo(Coordinate that) {

        int compare = this.start - that.start;
        
        if (compare == 0) {
            compare = this.stop - that.start;
        }
        
        return compare;
    }
}
