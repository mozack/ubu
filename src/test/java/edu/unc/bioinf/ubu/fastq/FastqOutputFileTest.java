package edu.unc.bioinf.ubu.fastq;

import static edu.unc.bioinf.ubu.fastq.FastqTestData.*;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.testng.annotations.Test;

/**
 * Unit tests for {@code FastqOutputFile}
 * 
 * @author lmose
 */
public class FastqOutputFileTest {
    
	private static final String OUTPUT_DIR = "src/test/java/edu/unc/bioinf/ubu/fastq/testdata/";
    private static final String OUTPUT_FILE = OUTPUT_DIR + "output.fastq";
    private static final String EXPECTED_OUTPUT_FILE = OUTPUT_DIR + "expected_output.fastq";

    @Test(groups = "unit")
    public void testOutputFile() throws IOException {
        File file = new File(OUTPUT_FILE);
        file.delete();
        assertFalse(file.exists(), "Precondition check failed.  output.fastq file could not be removed.");
        
        FastqOutputFile fastqFile = new FastqOutputFile();
        fastqFile.init(OUTPUT_FILE);
        fastqFile.write(REC1);
        fastqFile.write(REC2);
        fastqFile.close();

        assertEquals(loadFileContent(OUTPUT_FILE), loadFileContent(EXPECTED_OUTPUT_FILE));        
    }
    
    private String loadFileContent(String filename) throws IOException {
        
        StringBuffer content = new StringBuffer();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        
        String line = null;
        
        while ((line = reader.readLine()) != null) {
            content.append(line);
            content.append('\n');
        }
        
        return content.toString();
    }
}