package edu.unc.bioinf.ubu.fastq;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;
import static org.testng.Assert.fail;

import org.testng.annotations.Test;

/**
 * Unit tests for {@code FastqRecord}
 * 
 * @author lmose
 */
public class FastqRecordTest {
    
	private static final String FASTQ1_QUALITY_PHRED33 = "@@@DDBDDFFFFDEEE";
	private static final String FASTQ1_QUALITY_PHRED64 = "___ccacceeeecddd";
	
	
    private static final String[] FASTQ1 = new String[] {
        "@UNC16-SN851_55:2:1101:1487:1950/1",
        "CAACTAGATCAAAAAT",
        "+",
        FASTQ1_QUALITY_PHRED33
        };

    private static final String[] FASTQ2 = new String[] {
        "@UNC16-SN851_55:2:1101:1487:1950/2",
        "CGTATAGATCGGGGC",
        "+",
        "CGTATAGA"
        };
    
    private static final String[] FASTQ3 = new String[] {
        "@UNC16-SN851_55:2:1101:1487:2001/2",
        "CAACTAGATCAAAAAT",
        "+",
        "CCCFFFFEHH"
        };

    private static final String[] MISSING_SLASH_IN_ID = new String[] {
        "@UNC16-SN851_55:2:1101:1487:2001",
        "CAACTAGATCAAAAAT",
        "+",
        "CCCFFFFEHH"
        };
    
    private static final String[] FASTQ4 = new String[] {
        "@UNC15-SN850:105:D047RACXX:1:1101:1242:2131 1:N:0:ATCACG",
        "NACCTTTCTTGGCCAGGACTCGCAGTGGGTAACCTGTTTCATC",
        "+",
        "#1=DDFFFHHHHHJIJJIJGHHJJJGHJJDGIJJIJIIJJJJJ"
        };
    
    @Test(groups = "unit")
    public void testGetBaseId_slash() {
        FastqRecord rec = new FastqRecord(FASTQ1);
        assertEquals(rec.getBaseId(), "@UNC16-SN851_55:2:1101:1487:1950");
    }
    
    @Test(groups = "unit")
    public void testGetBaseId_space() {
        FastqRecord rec = new FastqRecord(FASTQ4);
        assertEquals(rec.getBaseId(), "@UNC15-SN850:105:D047RACXX:1:1101:1242:2131");
    }
    
    @Test(groups = "unit")
    public void testHasSameBaseId() {
        FastqRecord rec1 = new FastqRecord(FASTQ1);
        FastqRecord rec2 = new FastqRecord(FASTQ2);
        FastqRecord rec3 = new FastqRecord(FASTQ3);
        
        assertTrue(rec1.hasSameBaseId(rec2));
        assertFalse(rec2.hasSameBaseId(rec3));
    }
    
    @Test(groups = "unit")
    public void testInvalidNumLinesOnInput() {
        try {
            new FastqRecord(new String[] { "invalid1", "invalid2" });
            fail("Expected exception on invalid input to FastqRecord constructor");
        } catch (IllegalArgumentException e) {
            // OK - Exception is expected
        }
    }
    
    @Test(groups = "unit")
    public void testInvalidId() {
        FastqRecord rec = new FastqRecord(MISSING_SLASH_IN_ID);
        rec.getBaseId();
        assertEquals(rec.getId(), "@UNC16-SN851_55:2:1101:1487:2001");
    }
    
    @Test(groups = "unit")
    public void testStripNonReadInfoInId() {
        FastqRecord rec = new FastqRecord(FASTQ4.clone());
        rec.stripNonReadInfoInId();
        assertEquals(rec.getId(), "@UNC15-SN850:105:D047RACXX:1:1101:1242:2131"); 
    }
    
    @Test(groups = "unit")
    public void testAppendSuffixToId() {
        FastqRecord rec = new FastqRecord(FASTQ4.clone());
        rec.appendToId("/1");
        assertEquals(rec.getId(), "@UNC15-SN850:105:D047RACXX:1:1101:1242:2131 1:N:0:ATCACG/1");
    }
    
    @Test(groups = "unit")
    public void testPhred33To64() {
    	FastqRecord rec = new FastqRecord(FASTQ1.clone());
    	
    	assertEquals(FASTQ1_QUALITY_PHRED33, rec.getQuality());
    	rec.phred33To64();
    	assertEquals(FASTQ1_QUALITY_PHRED64, rec.getQuality());
    }
}