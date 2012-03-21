package edu.unc.bioinf.ubu.util;

import static junit.framework.Assert.assertEquals;

import org.testng.annotations.Test;

/**
 * Unit tests for {@code QualityConverter}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class QualityConverterTest {

	private static final String PHRED33 = "@@@DDBDDFFFFDEEE";
	private static final String PHRED64 = "___ccacceeeecddd";

	
    @Test (groups = "unit")
    public void testPhred33ToPhred64() {
    	QualityConverter converter = new QualityConverter();
    	assertEquals(PHRED64, converter.phred33ToPhred64(PHRED33));
    }

    @Test (groups = "unit")
    public void testPhred64ToPhred33() {
    	QualityConverter converter = new QualityConverter();
    	assertEquals(PHRED33, converter.phred64ToPhred33(PHRED64));
    }
}
