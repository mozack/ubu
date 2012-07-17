package edu.unc.bioinf.ubu.fastq;

import static org.testng.Assert.assertEquals;

import java.io.File;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.util.FileLoader;

/**
 * Test classes for {@code Sam2Fastq}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2FastqTest {
	
	private static final String TEST_DIR = "src/test/java/edu/unc/bioinf/ubu/fastq/testdata/";
	
	// Includes forward / reverse strand reads, unaligned reads and multi-mapper.
	private static final String TEST_PAIRED_END_SAM_FILE = TEST_DIR + "sam2fastq_test.sam";
	
	private static final String PAIRED_END_OUT1 = TEST_DIR + "sam2fastq_paired_out1.fastq";
	private static final String PAIRED_END_OUT2 = TEST_DIR + "sam2fastq_paired_out2.fastq";
	private static final String PAIRED_END_EXPECTED_OUT1 = TEST_DIR + "sam2fastq_expected_paired1.fastq";
	private static final String PAIRED_END_EXPECTED_OUT2 = TEST_DIR + "sam2fastq_expected_paired2.fastq";
	
	@BeforeMethod(groups = "unit")
	public void cleanup() {
		new File(PAIRED_END_OUT1).delete();
		new File(PAIRED_END_OUT2).delete();
	}
	
	@Test(groups = "unit")
	public void testSam2Fastq_pairedEndUsingReadIds() throws Exception {
		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.setEndSuffixes("/1", "/2");
		sam2Fastq.convert(TEST_PAIRED_END_SAM_FILE, PAIRED_END_OUT1, PAIRED_END_OUT2);
		
		FileLoader fileLoader = new FileLoader();
		assertEquals(fileLoader.loadFileContent(PAIRED_END_OUT1), fileLoader.loadFileContent(PAIRED_END_EXPECTED_OUT1));
		assertEquals(fileLoader.loadFileContent(PAIRED_END_OUT2), fileLoader.loadFileContent(PAIRED_END_EXPECTED_OUT2));
	}
}
