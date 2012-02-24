package edu.unc.bioinf.ubu.fastq;

import static edu.unc.bioinf.ubu.fastq.FastqTestData.*;
import static org.easymock.EasyMock.createMock;
import static org.easymock.EasyMock.expect;
import static org.easymock.EasyMock.replay;
import static org.easymock.EasyMock.verify;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * Unit tests for {@code FastqFilter}
 * 
 * @author lmose
 */
public class FastqFilterTest {

    private static final int MAX_CACHED_LINES = 3;
    
    private FastqInputFile input1;
    private FastqInputFile input2;
    private FastqOutputFile output1;
    private FastqOutputFile output2;
    private FastqFilter fastqFilter;
    private Object[] mocks;
    
    @SuppressWarnings("unused")
    @BeforeMethod(groups = "unit")
    private void setUp() throws FileNotFoundException, IOException {
        input1  = createMock(FastqInputFile.class);
        input2  = createMock(FastqInputFile.class);
        output1 = createMock(FastqOutputFile.class);
        output2 = createMock(FastqOutputFile.class);
        
        fastqFilter = new FastqFilter("input1", "input2", "output1", "output2", MAX_CACHED_LINES);
        fastqFilter.init(input1, input2, output1, output2);
        
        // Setup common expectations
        // "expect" is implied for void methods
        input1.init("input1", MAX_CACHED_LINES);
        input2.init("input2", MAX_CACHED_LINES);
        output1.init("output1");
        output2.init("output2");
        
        mocks = new Object[] { input1, input2, output1, output2 };
    }
    
    @Test(groups = "unit")
    public void testFilter_allReadsMatch() throws Exception {
        
        // Define input
        expect(input1.getRecord(1)).andReturn(REC1).anyTimes();
        expect(input1.getRecord(2)).andReturn(REC2).anyTimes();
        expect(input1.getRecord(3)).andReturn(REC3).anyTimes();
        expect(input1.getRecord(4)).andReturn(REC4).anyTimes();
        expect(input1.getRecord(5)).andReturn(REC5).anyTimes();
        expect(input1.getRecord(6)).andReturn(null).anyTimes();
        
        expect(input2.getRecord(1)).andReturn(REC1_2).anyTimes();
        expect(input2.getRecord(2)).andReturn(REC2_2).anyTimes();
        expect(input2.getRecord(3)).andReturn(REC3_2).anyTimes();
        expect(input2.getRecord(4)).andReturn(REC4_2).anyTimes();
        expect(input2.getRecord(5)).andReturn(REC5_2).anyTimes();
        expect(input2.getRecord(6)).andReturn(null).anyTimes();
        
        // Expected writes to output.
        output1.write(REC1);
        output1.write(REC2);
        output1.write(REC3);
        output1.write(REC4);
        output1.write(REC5);
        
        output2.write(REC1_2);
        output2.write(REC2_2);
        output2.write(REC3_2);
        output2.write(REC4_2);
        output2.write(REC5_2);
        
        output1.close();
        output2.close();
        
        replay(mocks);
        fastqFilter.filter();
        verify(mocks);
    }
    
    @Test(groups = "unit")
    public void testFilter_someReadsMatch() throws Exception {
        
        // Define input
        expect(input1.getRecord(1)).andReturn(REC1).anyTimes();
        expect(input1.getRecord(2)).andReturn(REC3).anyTimes();
        expect(input1.getRecord(3)).andReturn(REC4).anyTimes();
        expect(input1.getRecord(4)).andReturn(REC5).anyTimes();
        expect(input1.getRecord(5)).andReturn(null).anyTimes();
        
        expect(input2.getRecord(1)).andReturn(REC1_2).anyTimes();
        expect(input2.getRecord(2)).andReturn(REC2_2).anyTimes();
        expect(input2.getRecord(3)).andReturn(REC3_2).anyTimes();
        expect(input2.getRecord(4)).andReturn(REC5_2).anyTimes();
        expect(input2.getRecord(5)).andReturn(REC8_2).anyTimes();
        expect(input2.getRecord(6)).andReturn(null).anyTimes();
        
        // Expected writes to output.
        output1.write(REC1);
        output1.write(REC3);
        output1.write(REC5);
        
        output2.write(REC1_2);
        output2.write(REC3_2);
        output2.write(REC5_2);
        
        output1.close();
        output2.close();
        
        replay(mocks);
        fastqFilter.filter();
        verify(mocks);
    }
    
    @Test(groups = "unit")
    public void testFilter_noReadsMatch() throws Exception {
        
        // Define input
        expect(input1.getRecord(1)).andReturn(REC1).anyTimes();
        expect(input1.getRecord(2)).andReturn(REC2).anyTimes();
        expect(input1.getRecord(3)).andReturn(REC3).anyTimes();
        expect(input1.getRecord(4)).andReturn(REC4).anyTimes();
        expect(input1.getRecord(5)).andReturn(null).anyTimes();
        
        expect(input2.getRecord(1)).andReturn(REC5_2).anyTimes();
        expect(input2.getRecord(2)).andReturn(REC6_2).anyTimes();
        expect(input2.getRecord(3)).andReturn(REC7_2).anyTimes();
        expect(input2.getRecord(4)).andReturn(REC8_2).anyTimes();
        expect(input2.getRecord(5)).andReturn(null).anyTimes();

        // Output files should be closed with no writes.
        output1.close();
        output2.close();
        
        replay(mocks);
        fastqFilter.filter();
        verify(mocks);
    }
}