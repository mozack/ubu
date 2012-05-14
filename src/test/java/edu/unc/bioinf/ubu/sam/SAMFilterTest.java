package edu.unc.bioinf.ubu.sam;

import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * Unit tests for {@code SAMFilter}
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SAMFilterTest {
	
	private SAMFilter filter;
	
	@BeforeMethod(groups = "unit")
	@SuppressWarnings(value = "unused")
	private void setUp() {
		filter = new SAMFilter();
		filter.setPairedEnd(true);
		filter.setShouldStripIndels(true);
		filter.setMaxInsertLen(10000);
		filter.setMinMappingQuality(1);
	}
	
	private SAMRecord buildTestRead(String cigar, int mappingQuality, int insertLen) {
		SAMRecord read = new SAMRecord(new SAMFileHeader());
		
		read.setCigarString(cigar);
		read.setMappingQuality(mappingQuality);
		read.setInferredInsertSize(insertLen);
		
		return read;
	}

	@Test (groups = "unit")
	public void testIsReadIncluded_positive() {
		SAMRecord read = buildTestRead("50M", 1, 500);
		assertTrue(filter.isReadIncluded(read));
	}
	
	@Test (groups = "unit")
	public void testIsReadIncluded_excludeInsert() {
		SAMRecord read = buildTestRead("40M1I9M", 1, 500);
		assertFalse(filter.isReadIncluded(read));
	}
	
	@Test (groups = "unit")
	public void testIsReadIncluded_excludeLowMappingQuality() {
		SAMRecord read = buildTestRead("50M", 0, 500);
		assertFalse(filter.isReadIncluded(read));
	}

	@Test (groups = "unit")
	public void testIsReadIncluded_excludeLargeInsert() {
		SAMRecord read = buildTestRead("50M", 1, 100001);
		assertFalse(filter.isReadIncluded(read));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_positive() {
		SAMRecord read1 = buildTestRead("50M", 1, 500);
		SAMRecord read2 = buildTestRead("50M", 1, 500);
		assertTrue(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead1Delete() {
		SAMRecord read1 = buildTestRead("1D49M", 1, 500);
		SAMRecord read2 = buildTestRead("50M", 1, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead2Insert() {
		SAMRecord read1 = buildTestRead("50M", 1, 500);
		SAMRecord read2 = buildTestRead("49M1I", 1, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead1AndRead2Indels() {
		SAMRecord read1 = buildTestRead("10I10D30M", 1, 500);
		SAMRecord read2 = buildTestRead("30M10D10I", 1, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead1LowMappingQuality() {
		SAMRecord read1 = buildTestRead("50M", 0, 500);
		SAMRecord read2 = buildTestRead("50M", 1, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead2LowMappingQuality() {
		SAMRecord read1 = buildTestRead("50M", 1, 500);
		SAMRecord read2 = buildTestRead("50M", 0, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead1LargeInsert() {
		SAMRecord read1 = buildTestRead("50M", 1, 10001);
		SAMRecord read2 = buildTestRead("50M", 1, 500);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead2LargeInsert() {
		SAMRecord read1 = buildTestRead("50M", 1, 500);
		SAMRecord read2 = buildTestRead("50M", 1, 10001);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
	
	@Test (groups = "unit")
	public void testIsBasePairIncluded_excludeRead1LowMappingQualityRead2LargeInsert() {
		SAMRecord read1 = buildTestRead("50M", 0, 500);
		SAMRecord read2 = buildTestRead("50M", 1, 10001);
		assertFalse(filter.isReadPairIncluded(new ReadPair(read1, read2)));
	}
}
