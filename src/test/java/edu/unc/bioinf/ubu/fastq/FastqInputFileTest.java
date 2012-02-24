package edu.unc.bioinf.ubu.fastq;

import static edu.unc.bioinf.ubu.fastq.FastqTestData.*;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.fail;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * Unit tests for {@code FastqInputFile}
 * 
 * @author lmose
 */
public class FastqInputFileTest {
    
    private static final String TEST_FASTQ_FILE = "src/test/java/edu/unc/bioinf/ubu/fastq/testdata/input.fastq";
    private static final int MAX_LINES_TO_CACHE = 3;
    
    private FastqInputFile file;
    
    @BeforeMethod(groups = "unit")
    @SuppressWarnings("unused")
    private void setUp() throws FileNotFoundException {
        file = new FastqInputFile();
        file.init(TEST_FASTQ_FILE, MAX_LINES_TO_CACHE);
    }

    @AfterMethod(groups = "unit")
    @SuppressWarnings("unused")
    private void tearDown() throws IOException {
        if (file != null) {
            file.close();
        }
    }

    @Test(groups = "unit")
    public void testGetRecord() throws Exception {
        FastqRecord rec;
        
        
        // Read first record 
        rec = file.getRecord(1);
        assertEquals(rec, REC1);
        
        rec = file.getRecord(2);
        assertEquals(rec, REC2);
        
        rec = file.getRecord(3);
        assertEquals(rec, REC3);

        // Go back a record
        rec = file.getRecord(2);
        assertEquals(rec, REC2);

        // Re-read
        rec = file.getRecord(3);
        assertEquals(rec, REC3);
        
        // Go forward
        rec = file.getRecord(4);
        assertEquals(rec, REC4);
        
        // Last record
        rec = file.getRecord(5);
        assertEquals(rec, REC5);
        
        // End of file
        rec = file.getRecord(6);
        assertEquals(rec, null);
    }
    
    @Test(groups = "unit")
    public void testGetRecord_indexLessThanOne() throws Exception {
        
        try {
            file.getRecord(0);
            fail("Expected exception on attempt read prior to the first record.");
        } catch (IllegalArgumentException e) {
            // OK - Exception is expected
        }
    }
    
    @Test(groups = "unit")
    public void testGetRecord_readBackBeforeCache() throws Exception {
        try {
            file.getRecord(5);
            file.getRecord(1);
            
            fail("Expected exception on attempt to 'read back' greater than cached record amount");
        } catch (UnsupportedOperationException e) {
            // OK - Exception is expected
        }
    }
    
    @Test(groups = "unit")
    public void testGetNextRecord() throws Exception {
        
        if (file != null) {
            file.close();
        }
        
        file = new FastqInputFile();
        file.init(TEST_FASTQ_FILE);
        
        FastqRecord rec;
        
        rec = file.getNextRecord();
        assertEquals(rec, REC1);

        rec = file.getNextRecord();
        assertEquals(rec, REC2);
        
        rec = file.getNextRecord();
        assertEquals(rec, REC3);
        
        rec = file.getNextRecord();
        assertEquals(rec, REC4);
        
        rec = file.getNextRecord();
        assertEquals(rec, REC5);
        
        rec = file.getNextRecord();
        assertEquals(rec, null);
    }
}
