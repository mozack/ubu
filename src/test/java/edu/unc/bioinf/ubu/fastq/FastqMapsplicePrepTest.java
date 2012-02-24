package edu.unc.bioinf.ubu.fastq;

import static org.easymock.EasyMock.expect;
import static org.easymock.EasyMock.replay;
import static org.easymock.EasyMock.verify;

import org.easymock.EasyMock;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class FastqMapsplicePrepTest {
    
    private FastqInputFile input;
    private FastqOutputFile output;
    private FastqMapsplicePrep mapsplicePrep;
    
    @BeforeMethod(groups = "unit")
    @SuppressWarnings("unused")
    private void setUp() {
        input = EasyMock.createMock(FastqInputFile.class);
        output = EasyMock.createMock(FastqOutputFile.class);
        mapsplicePrep = new FastqMapsplicePrep(input, output, "/1");
    }

    @Test(groups = "unit")
    public void testProcess_casava18() throws Exception {
       
        // Setup expectations
        expect(input.getNextRecord()).andReturn(FastqTestData.CASAVA_1_8_REC1);
        output.write(FastqTestData.CASAVA_CONVERTED_REC1);
        expect(input.getNextRecord()).andReturn(FastqTestData.CASAVA_1_8_REC2);
        output.write(FastqTestData.CASAVA_CONVERTED_REC2);
        expect(input.getNextRecord()).andReturn(FastqTestData.CASAVA_1_8_REC3);
        output.write(FastqTestData.CASAVA_CONVERTED_REC3);
        expect(input.getNextRecord()).andReturn(null);
        input.close();
        output.close();

        // Execute and verify output
        replay(input, output);
        mapsplicePrep.process();
        verify(input, output);
    }
    
    @Test(groups = "unit")
    public void testProcess_preCasava18() throws Exception {
        
    }
}
