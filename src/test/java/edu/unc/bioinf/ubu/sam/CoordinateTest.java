package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.Coordinate;

/**
 * Unit tests for {@code Coordinate}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class CoordinateTest {
    
    private Coordinate c = new Coordinate(100, 200);

    @Test (groups = "unit")
    public void testContains() {
        
        assertFalse(c.contains(99));
        assertTrue(c.contains(100));
        assertTrue(c.contains(101));
        assertTrue(c.contains(150));
        assertTrue(c.contains(199));
        assertTrue(c.contains(200));
        assertFalse(c.contains(201));
    }
    
    @Test (groups = "unit")
    public void testGetLength() {
        // Coordinate Start and stop are inclusive
        // 100 -> 200 = 101
        assertEquals(c.getLength(), 101);
    }
}
