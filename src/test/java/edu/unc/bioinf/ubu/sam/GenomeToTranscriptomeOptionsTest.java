package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.GenomeToTranscriptomeOptions;

/**
 * Unit tests for {@code GenomeToTranscriptomeOptions}
 * @author lmose
 */
public class GenomeToTranscriptomeOptionsTest {

    private GenomeToTranscriptomeOptions parser;
    
    @BeforeMethod (groups = "unit")
    void setUp() {
        parser = new GenomeToTranscriptomeOptions();
    }
    
    @Test (groups = "unit")
    public void testNoParams() {
        parser.parseOptions(
                "".split(" "));
            
        assertFalse(parser.isValid());
    }
    
    @Test (groups = "unit")
    public void testHelp() {
        parser.parseOptions(
            new String[] {
                "--help" 
            });
        
        assertFalse(parser.isValid());
    }

    @Test (groups = "unit")
    public void testMinParams() {
        parser.parseOptions(
                "--bed bedfile --in infile --out outfile".split(" "));
            
        assertTrue(parser.isValid());
        validateBasicParams();
        assertTrue(parser.isPositiveStrandReportingOnly());
        assertFalse(parser.shouldOutputXgTags());
    }
    
    @Test (groups = "unit")
    public void testAllParams() {
        parser.parseOptions(
                "--bed bedfile --in infile --out outfile --offset 50 --dupes dupefile --order orderfile --reverse --xgtags".split(" "));
            
        assertTrue(parser.isValid());
        validateBasicParams();
        validateDupes();
        validateOrderFile();
        assertEquals(parser.getReadOffset(), 50);
        assertFalse(parser.isPositiveStrandReportingOnly());
        assertTrue(parser.shouldOutputXgTags());
    }
    
    @Test (groups = "unit")
    public void testInitialMapspliceParams() {
        parser.parseOptions(
                "--bed bedfile --in infile --out outfile --dupes dupefile --order orderfile --xgtags".split(" "));

        assertTrue(parser.isValid());
        validateBasicParams();
        validateDupes();
        validateOrderFile();
        assertEquals(parser.getReadOffset(), 25);
        assertTrue(parser.isPositiveStrandReportingOnly());
        assertTrue(parser.shouldOutputXgTags());
    }
    
    @Test (groups = "unit")
    public void testMissingFile() {
        parser.parseOptions(
                "--bed bedfile --in infile --out ".split(" "));
        
        assertFalse(parser.isValid());
    }

    private void validateBasicParams() {
        assertEquals(parser.getBedFile(), "bedfile");
        assertEquals(parser.getInputAlignmentFile(), "infile");
        assertEquals(parser.getOutputAlignmentFile(), "outfile");
    }
    
    private void validateDupes() {
        assertEquals(parser.getDuplicatesFile(), "dupefile");
    }
    
    private void validateOrderFile() {
        assertEquals(parser.getOrderingFastaFile(), "orderfile");
    }
}
