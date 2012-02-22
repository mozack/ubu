package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertEquals;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.ReverseComplementor;

/**
 * Unit tests for {@code ReverseComplementor}
 */
public class ReverseComplementorTest {
    
    private byte[] BASES;
    
    @BeforeMethod (groups = "unit")
    public void setUp() {
        // Re-initialize bases each time in case of mutations.
        BASES = "CTCTNCTGATCCCCACCTCCAAATATCTCATCAACAACCAACTAATCACC".getBytes();
    }

    @Test (groups = "unit")
    public void testReverse() {
        ReverseComplementor rc = new ReverseComplementor();
        byte[] reverse = rc.reverse(BASES);
        
        assertEquals("CCACTAATCAACCAACAACTACTCTATAAACCTCCACCCCTAGTCNTCTC".getBytes(), reverse);
        // input should be unchanged
        assertEquals("CTCTNCTGATCCCCACCTCCAAATATCTCATCAACAACCAACTAATCACC".getBytes(), BASES);
    }
    
    @Test (groups = "unit")
    public void testComplement() {
        ReverseComplementor rc = new ReverseComplementor();
        rc.complementInPlace(BASES);
        
        assertEquals("GAGANGACTAGGGGTGGAGGTTTATAGAGTAGTTGTTGGTTGATTAGTGG".getBytes(), BASES);
    }
    
    @Test (groups = "unit")
    public void testReverseComplement() {
        ReverseComplementor rc = new ReverseComplementor();
        byte[] reverseComplement = rc.reverseComplement(BASES);
        
        assertEquals("GGTGATTAGTTGGTTGTTGATGAGATATTTGGAGGTGGGGATCAGNAGAG".getBytes(), reverseComplement);
        
        // input should be unchanged
        assertEquals("CTCTNCTGATCCCCACCTCCAAATATCTCATCAACAACCAACTAATCACC".getBytes(), BASES);
        
    }
}
