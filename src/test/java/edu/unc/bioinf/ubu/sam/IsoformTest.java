package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.Coordinate;
import edu.unc.bioinf.ubu.sam.Isoform;

/**
 * Unit tests for {@code Isoform}
 * 
 * @author lmose
 */
public class IsoformTest {
    
    private static final List<Coordinate> exons1 = Arrays.asList(
                new Coordinate(1001, 1300),
                new Coordinate(1501, 1700),
                new Coordinate(1901, 2000));
    
    private static final Isoform iso1 = new Isoform(
            "iso1", new Coordinate(1001,2000), "+", exons1);

    private static final List<Coordinate> exons2 = Arrays.asList(
            new Coordinate(10000, 10999));
    
    private static final Isoform iso2 = new Isoform(
            "iso1", new Coordinate(10000, 10999), "-", exons2);


    @Test (groups = "unit")
    public void testContainsWithinGenomicRange() {
        assertTrue(iso1.containsWithinGenomicRange(1001, 2000));
        assertTrue(iso1.containsWithinGenomicRange(1400, 1600));
        
        assertFalse(iso1.containsWithinGenomicRange(900, 1000));
        assertFalse(iso1.containsWithinGenomicRange(1000, 1001));
        assertFalse(iso1.containsWithinGenomicRange(1000, 1100));
        assertFalse(iso1.containsWithinGenomicRange(1000, 2001));
        assertFalse(iso1.containsWithinGenomicRange(2000, 2001));
        assertFalse(iso1.containsWithinGenomicRange(2001, 2100));
    }
    
    @Test (groups = "unit")
    public void testGetLength() {
        // [1001,1300] = 300
        // [1501,1700] = 200
        // [1901,2000] = 100
        assertEquals(iso1.getLength(), 600);
    }
    
    @Test (groups = "unit")
    public void testIsPositiveStrand() {
        assertTrue(iso1.isPositiveStrand());
        assertFalse(iso2.isPositiveStrand());
    }
    
    @Test (groups = "unit")
    public void testIsNegativeStrand() {
        assertFalse(iso1.isNegativeStrand());
        assertTrue(iso2.isNegativeStrand());
    }
    
    @Test (groups = "unit")
    public void testMatch_startOfIsoform() {
        SAMRecord read = createRead("chr1", 1001, "50M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 1);
        assertEquals(coords.get(0).getStart(), 1);
        assertEquals(coords.get(0).getStop(), 50);
    }
    
    @Test (groups = "unit")
    public void testMatch_endOfIsoform() {
        SAMRecord read = createRead("chr1", 1951, "50M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 1);
        assertEquals(coords.get(0).getStart(), 551);
        assertEquals(coords.get(0).getStop(), 600);
    }
    
    @Test (groups = "unit")
    public void testMatch_middleExon() {
        SAMRecord read = createRead("chr1", 1600, "50M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 1);
        assertEquals(coords.get(0).getStart(), 400);
        assertEquals(coords.get(0).getStop(), 449);
    }
    
    @Test (groups = "unit")
    public void testMatch_readSkipsFirstIntron() {        
        SAMRecord read = createRead("chr1", 1291, "10M200N40M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 291);
        assertEquals(coords.get(0).getStop(), 300);
        assertEquals(coords.get(1).getStart(), 301);
        assertEquals(coords.get(1).getStop(), 340);
    }
    
    @Test (groups = "unit")
    public void testMatch_readSkipsSecondIntron() {
        SAMRecord read = createRead("chr1", 1700, "1M200N49M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 500);
        assertEquals(coords.get(0).getStop(), 500);
        assertEquals(coords.get(1).getStart(), 501);
        assertEquals(coords.get(1).getStop(), 549);
    }
    
    @Test (groups = "unit")
    public void testMatch_readSpansMultipleIntrons() {
        SAMRecord read = createRead("chr1", 1291, "10M200N200M200N15M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 3);
        assertEquals(coords.get(0).getStart(), 291);
        assertEquals(coords.get(0).getStop(), 300);
        assertEquals(coords.get(1).getStart(), 301);
        assertEquals(coords.get(1).getStop(), 500);
        assertEquals(coords.get(2).getStart(), 501);
        assertEquals(coords.get(2).getStop(), 515);
    }
    
    @Test (groups = "unit")
    public void testMatch_insert() {
        SAMRecord read = createRead("chr1", 1001, "33M1I16M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 1);
        assertEquals(coords.get(0).getStop(), 33);
        assertEquals(coords.get(1).getStart(), 34);
        assertEquals(coords.get(1).getStop(), 49);
    }
    
    @Test (groups = "unit")
    public void testMatch_insertReadAtEndEdgeOfExon() {
        SAMRecord read = createRead("chr1", 1652, "33M1I16M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 452);
        assertEquals(coords.get(0).getStop(), 484);
        assertEquals(coords.get(1).getStart(), 485);
        assertEquals(coords.get(1).getStop(), 500);
    }
    
    @Test (groups = "unit")
    public void testMatch_delete() {
        SAMRecord read = createRead("chr1", 1001, "49M2D1M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 1);
        assertEquals(coords.get(0).getStop(), 49);
        assertEquals(coords.get(1).getStart(), 52);
        assertEquals(coords.get(1).getStop(), 52);
    }
    
    @Test (groups = "unit")
    public void testMatch_deleteMovesReadIntoIntron() {        
        SAMRecord read = createRead("chr1", 1250, "49M2D1M");
        
        List<Coordinate> coords = iso1.match(read);
        assertTrue(coords.isEmpty());
    }
    
    @Test (groups = "unit")
    public void testMatch_deleteReadAtEndEdgeOfExon() {
        SAMRecord read = createRead("chr1", 1249, "49M2D1M");
        
        List<Coordinate> coords = iso1.match(read);
        assertEquals(coords.size(), 2);
        assertEquals(coords.get(0).getStart(), 249);
        assertEquals(coords.get(0).getStop(), 297);
        assertEquals(coords.get(1).getStart(), 300);
        assertEquals(coords.get(1).getStop(), 300);
    }
    
    @Test (groups = "unit")
    public void testMatch_readInStartOfFirstIntron() {
        
        SAMRecord read = createRead("chr1", 1291, "11M199N40M");
        
        List<Coordinate> coords = iso1.match(read);
        assertTrue(coords.isEmpty());
    }
    
    @Test (groups = "unit")
    public void testMatch_readInEndOfFirstIntron() {
        
        SAMRecord read = createRead("chr1", 1291, "10M201N40M");
        
        List<Coordinate> coords = iso1.match(read);
        assertTrue(coords.isEmpty());
    }
    
    private SAMRecord createRead(String refName, int alignmentStart, String cigar) {
        SAMRecord read = new SAMRecord(new SAMFileHeader());
        
        read.setReferenceName(refName);
        read.setAlignmentStart(alignmentStart);
        read.setCigarString(cigar);

        return read;
    }
}
