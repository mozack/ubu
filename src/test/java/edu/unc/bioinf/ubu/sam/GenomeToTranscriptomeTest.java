package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * System test for {@code GenomeToTranscriptome}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class GenomeToTranscriptomeTest {
	
	private static final String TESTDATA_DIR = "src/test/java/edu/unc/bioinf/ubu/sam/testdata/";
	private static final String OUTPUT_DIR = TESTDATA_DIR + "test_output/";
	private static final String BED = TESTDATA_DIR + "unc_hg19.bed";
	private static final String TRANSCRIPTS = TESTDATA_DIR + "hg19_M_rCRS_ref.transcripts.fa";
	private static final String INPUT = TESTDATA_DIR + "system_test_sorted_by_ref_and_read.bam";
	private static final String OUTPUT = OUTPUT_DIR + "GenomeToTranscriptomeTest.bam";
	private static final String DUPES = OUTPUT_DIR + "dupes.txt";
	private static final String EXPECTED_OUTPUT = TESTDATA_DIR + "transcriptome.bam";

    @BeforeMethod (groups = "system")
    void setUp() throws IOException, InterruptedException {
    	File output = new File(OUTPUT);
    	output.delete();
    	
    	File dupes = new File(DUPES);
    	dupes.delete();
    	
    	String rmTestFiles = "rm " + OUTPUT_DIR + "*";
    	System.out.println(rmTestFiles);
    	Process p = Runtime.getRuntime().exec(rmTestFiles);
    	p.waitFor();
    	
        File outputDir = new File(OUTPUT_DIR);
        outputDir.delete();
        
        boolean isDirCreated = outputDir.mkdir();
        if (!isDirCreated) {
        	throw new RuntimeException("Unable to create directory: " + OUTPUT_DIR);
        }
    }
    
    @Test (groups = "system")
    public void testGenomeToTranscriptome() throws Exception {
        String argz = 
            "--bed " + BED +
            " --offset 300 " +
            " --order " + TRANSCRIPTS +
            " --xgtags " +
            " --reverse " +
            " --in " + INPUT +
            " --out " + OUTPUT
            ;
        
        GenomeToTranscriptome.run(argz.split(" "));
        
        String command = "diff " + OUTPUT + " " + EXPECTED_OUTPUT;
        System.out.println("Executing: " + command);
        Process p = Runtime.getRuntime().exec(command);
        
        int diff = p.waitFor();
        if (diff != 0) {
        	Assert.fail("Expected transcriptome output does not match.  diff returned: " + diff);
        }
    }
}
