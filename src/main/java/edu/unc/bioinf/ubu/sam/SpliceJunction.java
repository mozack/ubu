package edu.unc.bioinf.ubu.sam;

/**
 * Represents a splice junction
 * TODO: Rename this to something more generic.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SpliceJunction implements Comparable<SpliceJunction> {

    private final String chromosome;
    private final Coordinate coordinate;
    
    public SpliceJunction(String chromosome, int start, int stop) {
        this.chromosome = chromosome;
        this.coordinate = new Coordinate(start, stop);
    }
    
    @Override
    public boolean equals(Object obj) {
        SpliceJunction that = (SpliceJunction) obj;
        
        return this.chromosome.equals(that.chromosome) && this.coordinate.equals(that.coordinate);
    }
    
    @Override
    public int hashCode() {
        int result = 17;
        
        result = 31 * result + chromosome.hashCode();
        result = 31 * result + coordinate.hashCode();
        
        return result;
    }

    @Override
    public int compareTo(SpliceJunction that) {
        int compare = this.chromosome.compareTo(that.chromosome);
        
        if (compare == 0) {
            compare = this.coordinate.compareTo(that.coordinate);
        }
        
        return compare;
    }
    
    public String toString() {
        return chromosome + ":" + coordinate;
    }
}
