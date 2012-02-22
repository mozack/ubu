package edu.unc.bioinf.ubu.sam;

import java.io.File;

import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.BedReader;

/**
 * 
 * @author lmose
 */
public class GenomeToTranscriptomeConverterTest {

    private static final String BED_FILE = "java/test/net/sourceforge/seqware/pipeline/util/ucsc_known_gene_bed.txt";
    private static final String INPUT1_BAM = "java/test/net/sourceforge/seqware/pipeline/util/input1.bam";
    
    private static final String OUTPUT1_BAM = "java/test/net/sourceforge/seqware/pipeline/util/output1.bam";
    private static final String OUTPUT1_SAM = "java/test/net/sourceforge/seqware/pipeline/util/output1.sam";
    private static final String OUTPUT1_SAM_HEADER = "java/test/net/sourceforge/seqware/pipeline/util/output1_header.sam";
    private static final String OUTPUT1_SORTED = "java/test/net/sourceforge/seqware/pipeline/util/output1_sorted";
    
    private static final String[] OUTPUT_FILES = new String[] {
        OUTPUT1_BAM,
        OUTPUT1_SAM,
        OUTPUT1_SORTED
    };
    
    @BeforeMethod(groups = "unit")
    @AfterMethod(groups = "unit")
    public void cleanUp() throws Exception {
        removeFiles(OUTPUT_FILES);
    }
    
    public void removeFiles(String[] files) {
        for (String filename : files) {
            File file = new File(filename);
            file.delete();
        }
    }
    
    // TODO: This isn't really a unit test
//    @Test(groups = "unit")
    public void testConvert() throws Exception {
        BedReader bedReader = new BedReader();
        bedReader.buildReadToIsoformIndex(BED_FILE, 200);
        
//        GenomeToTranscriptomeConverter converter = new GenomeToTranscriptomeConverter(bedReader);
//        
//        converter.convertFile(INPUT1_BAM, OUTPUT1_BAM);
        
        // Sort and convert the output bam to a SAM
        Runtime.getRuntime().exec("samtools sort " + OUTPUT1_BAM + " " + OUTPUT1_SORTED);
        Runtime.getRuntime().exec("samtools view " + OUTPUT1_SORTED + " > " + OUTPUT1_SAM);
        Runtime.getRuntime().exec("samtools view -H " + OUTPUT1_SORTED + " > " + OUTPUT1_SAM_HEADER);
        
        
    }
}
