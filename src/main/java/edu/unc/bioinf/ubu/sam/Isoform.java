package edu.unc.bioinf.ubu.sam;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Representation of an Isoform.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Isoform {
    
    public enum Strand {
        FORWARD,
        REVERSE
    };
    
    private static final String FORWARD_STRAND = "+";
    private static final String REVERSE_STRAND = "-";

    private String isoformId;
    private Coordinate genomicCoord;
    private List<Coordinate> exons;
    private Strand strand;
    private int length;
    
    public Isoform(String isoformId, Coordinate genomicCoord, String strand, List<Coordinate> exons) {
        this.isoformId = isoformId;
        this.genomicCoord = genomicCoord;
        this.exons = exons;
        this.length = calcLength(); 
        
        if (strand.equals(FORWARD_STRAND)) {
            this.strand = Strand.FORWARD;
        } else if (strand.equals(REVERSE_STRAND)) {
            this.strand = Strand.REVERSE;
        } else {
            throw new IllegalArgumentException("Invalid strand: " + strand + " for isoform: " + isoformId);
        }
    }

    /**
     * Returns true if the input coordinates are within this isoform's
     * genomic coordinates. 
     */
    public boolean containsWithinGenomicRange(int genomicStartPos, int genomicEndPos) {
        return genomicCoord.contains(genomicStartPos) && genomicCoord.contains(genomicEndPos);
    }
    
    public boolean isNegativeStrand() {
        return strand == Strand.REVERSE;
    }
    
    public boolean isPositiveStrand() {
        return strand == Strand.FORWARD;
    }
    
    public int getLength() {
        return length;
    }
    
    private int calcLength() {
        int length = 0;
        
        for (Coordinate exon : exons) {
            length += exon.getLength();
        }
        
        return length;
    }
    
    public Coordinate getGenomicCoordinate() {
        return this.genomicCoord;
    }
    
    public String getIsoformId() {
        return isoformId;
    }
    
    public String toString() {
        return isoformId + "," + genomicCoord;
    }
    
    private Coordinate getTranscriptCoordinates(int blockIsoformStart, int blockIsoformStop) {
        return new Coordinate(blockIsoformStart, blockIsoformStop);
    }
    

    /**
     * Attempts to match the input read with this Isoform.
     * If it is determined that the read is a match, a list of 1 based isoform coordinates
     * is returned.  A separate coordinate is generated for each alignment block in the read.
     * If the read is not a match, an empty list is returned.
     * <br>
     * ASSUMES CHROMOSOME HAS ALREADY BEEN MATCHED TO THIS ISOFORM!
     */
    public List<Coordinate> match(SAMRecord read) {
        
        List<Coordinate> transcriptCoordinates = new ArrayList<Coordinate>();
        
        boolean isFirstBlock = true;
        int prevBlockIsoformStop = 0;
        
        // For each block
        for (ReadBlock block : ReadBlock.getReadBlocks(read)) {
            
            if (block.getType() == CigarOperator.M) {
                boolean isBlockMatch = false;
                
                int blockGenomeStart = block.getReferenceStart();
                // i.e. start = 1, len = 50, so stop = 1 + 50 - 1 = 50
                int blockGenomeStop  = blockGenomeStart + block.getLength() - 1;
                
                int exonStartInIsoform = 1;
                
                // Find matching exon for this block
                for (Coordinate exon : exons) {
                    if ((exon.contains(blockGenomeStart)) && (exon.contains(blockGenomeStop))) {
                        // (1 based coordinate) - (1 based coordinate) + (offset into isoform)
                        int blockIsoformStart = blockGenomeStart - exon.getStart() + exonStartInIsoform ;
                        int blockIsoformStop  = blockGenomeStop - exon.getStart() + exonStartInIsoform;
                        
                        // If this is the first block, we have a match.
                        // For non-first blocks, if this block starts immediately after the previous
                        // block's stop (, we have a match
                        if ((isFirstBlock) || (blockIsoformStart == prevBlockIsoformStop + 1)) {
                            transcriptCoordinates.add(getTranscriptCoordinates(blockIsoformStart, blockIsoformStop));
                            isBlockMatch = true;
                            prevBlockIsoformStop = blockIsoformStop;
                            break;
                        }
                    }
                    
                    exonStartInIsoform += exon.getLength();
                }
                
                // If any block doesn't match, this read is not a match
                if (isBlockMatch == false) {
                    return Collections.EMPTY_LIST;
                }
                
                isFirstBlock = false;
            } else if (block.getType() == CigarOperator.D) {
                prevBlockIsoformStop += block.getLength();
            } else if (block.getType() == CigarOperator.I) {
                // Insert requires no change.
            } else if (block.getType() == CigarOperator.S) {
                // This appears to indicate a fusion for Mapsplice.  This won't map to a transcript,
                // so discard it.
                return Collections.EMPTY_LIST;
            }
        }
        
        return transcriptCoordinates;
    }

    /**
     * Simple Comparator class that can be used to sort Isoforms by the order
     * defined by the IsoformOrderLoader
     */
    public static class IsoformOrderComparator implements Comparator<Isoform> {
        
        private IsoformOrderLoader isoformOrderLoader;

        public IsoformOrderComparator(IsoformOrderLoader isoformOrderLoader) {
            this.isoformOrderLoader = isoformOrderLoader;
        }
        
        @Override
        public int compare(Isoform isoform1, Isoform isoform2) {
            return isoformOrderLoader.getOrder(isoform1.getIsoformId()).compareTo(
                   isoformOrderLoader.getOrder(isoform2.getIsoformId()));
        }
    }
}
