package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMRecord;

import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.sam.ReadPair;
import edu.unc.bioinf.ubu.sam.SamReadPairReader;

/**
 * Unit tests for {@code SamReadPairReader}
 * 
 * @author lmose
 */
public class SamReadPairReaderTest {

    private static final String TEST1_SAM_FILE = "src/test/java/edu/unc/bioinf/ubu/sam/testdata/test1.sam";
    private static final String MULTI_MATE_SAM_FILE = "src/test/java/edu/unc/bioinf/ubu/sam/testdata/multi_mate.sam";
    
    @Test (groups = "unit")
    public void testReadBasicSam() {
        SamReadPairReader reader = new SamReadPairReader(TEST1_SAM_FILE);
        
        List<ReadPair> readPairs = new ArrayList<ReadPair>();
        
        for (ReadPair readPair : reader) {
            readPairs.add(readPair);
        }
        
        assertEquals(3, readPairs.size());
        //UNC11-SN627_70:1:1101:1075:187685/1     99      chrM    8610    33      50M     =       8671
        verifyRead(readPairs.get(0).getRead1(), "UNC11-SN627_70:1:1101:1075:187685/1", "chrM", 8610, "50M", "chrM", 8671);
        //UNC11-SN627_70:1:1101:1075:187685/2     147     chrM    8671    55      50M     =       8610
        verifyRead(readPairs.get(0).getRead2(), "UNC11-SN627_70:1:1101:1075:187685/2", "chrM", 8671, "50M", "chrM", 8610);
        
        //UNC11-SN627_70:1:1101:4231:57374/1  99  chr1    144829336   60  50M =   144829414   127 CTTCAGCCTCCAATT
        verifyRead(readPairs.get(1).getRead1(), "UNC11-SN627_70:1:1101:4231:57374/1", "chr1", 144829336, "50M", "chr1", 144829414);
        //UNC11-SN627_70:1:1101:4231:57374/2  147 chr1    144829414   0   50M =   144829336   -127    CTCCTCTCAATTCCA
        verifyRead(readPairs.get(1).getRead2(), "UNC11-SN627_70:1:1101:4231:57374/2", "chr1", 144829414, "50M", "chr1", 144829336);

        // Read 1 is always the read with the first position.
        // So, the following appears reversed.
        
        //UNC11-SN627_70:1:1101:4231:57374/2    417 chr1    146033486   0   50M =   146033564   127 CTGCCAGCGGCTCCT
        verifyRead(readPairs.get(2).getRead1(), "UNC11-SN627_70:1:1101:4231:57374/2", "chr1", 146033486, "50M", "chr1", 146033564);
        //UNC11-SN627_70:1:1101:4231:57374/1  339 chr1    146033564   60  50M =   146033486   -127    TCC1TCTAGCTGTGG
        verifyRead(readPairs.get(2).getRead2(), "UNC11-SN627_70:1:1101:4231:57374/1", "chr1", 146033564, "50M", "chr1", 146033486);
    }
    
    @Test (groups = "unit")
    public void testMultipleMates() {
        SamReadPairReader reader = new SamReadPairReader(MULTI_MATE_SAM_FILE);
        
        List<ReadPair> readPairs = new ArrayList<ReadPair>();
        
        for (ReadPair readPair : reader) {
            readPairs.add(readPair);
        }
        
        assertEquals(readPairs.size(), 2);
        
        // UNC12-SN629_146:2:1101:2003:111782/1	323	chr22	39712776	60	50M	=	39712828	-3	GCTGCTTCTTGCCATCCTCATCCTGCCATTTCTTGCAGTACTTGGTAAAG	CCCFFFFFHHHHHJJJJIJJJJJJJJJIJJJJJJJJJJJJJJJJJFHIIJ	NM:i:0	IH:i:3	HI:i:2
        verifyRead(readPairs.get(0).getRead1(), "UNC12-SN629_146:2:1101:2003:111782/1", "chr22", 39712776, "50M", "chr22", 39712828);
        // UNC12-SN629_146:2:1101:2003:111782/2	163	chr22	39712828	57	19M619N31M	chr6	31248471	8465025	CTTCTTCTTAGATTTATGCCAATTCTTATAGAAACGCCTCTTGCATTCAA	CCCFFFFFHHHHHIJJJJJJJJJJIIJJJJJJJIJIIJIJJJJJGEIJIJ	NM:i:1	XS:A:-	XF:Z:CTAC,	IH:i:1	HI:i:1
        verifyRead(readPairs.get(0).getRead2(), "UNC12-SN629_146:2:1101:2003:111782/2", "chr22", 39712828, "19M619N31M", "chr6", 31248471);
        
        // UNC12-SN629_146:2:1101:2003:111782/1	339	chr6	31248471	60	50M	chr22	39712828	-8464308	CTTTACCAAGTACTGCAAGAAATGGCAGGATGAGGATGGCAAGAAGCAGC	JIIHFJJJJJJJJJJJJJJJJJIJJJJJJJJJIJJJJHHHHHFFFFFCCC	NM:i:0	IH:i:3	HI:i:3        
        verifyRead(readPairs.get(1).getRead1(), "UNC12-SN629_146:2:1101:2003:111782/1", "chr6", 31248471, "50M", "chr22", 39712828);
        // UNC12-SN629_146:2:1101:2003:111782/2	163	chr22	39712828	57	19M619N31M	chr6	31248471	8465025	CTTCTTCTTAGATTTATGCCAATTCTTATAGAAACGCCTCTTGCATTCAA	CCCFFFFFHHHHHIJJJJJJJJJJIIJJJJJJJIJIIJIJJJJJGEIJIJ	NM:i:1	XS:A:-	XF:Z:CTAC,	IH:i:1	HI:i:1
        verifyRead(readPairs.get(1).getRead2(), "UNC12-SN629_146:2:1101:2003:111782/2", "chr22", 39712828, "19M619N31M", "chr6", 31248471);        
    }

    private void verifyRead(SAMRecord read, String name, String ref, int position, String cigar, String mateRef, int matePosition) {
        assertEquals(name, read.getReadName());
        assertEquals(ref, read.getReferenceName());
        assertEquals(position, read.getAlignmentStart());
        assertEquals(cigar, read.getCigarString());
        assertEquals(mateRef, read.getMateReferenceName());
        assertEquals(matePosition, read.getMateAlignmentStart());
    }
}
