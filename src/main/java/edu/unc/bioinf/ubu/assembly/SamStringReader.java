package edu.unc.bioinf.ubu.assembly;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class SamStringReader {

	public SAMRecord getRead(String record) {
		InputStream is = new ByteArrayInputStream(record.getBytes());
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader reader = new SAMFileReader(is);
		
		SAMRecord read = reader.iterator().next();
		
		try {
			is.close();
		} catch (IOException e) {
			// This should not happen, so convert to RuntimeException
			e.printStackTrace();
			throw new RuntimeException(e);
		} finally {
			reader.close();
		}
		
		return read;
	}
	
	/*
	public static void main(String[] args) {
//		String readStr = "UNC11-SN627:184:81RL6ABXX:8:2103:10602:9739:2	16	chr11	118380105	37	100M	*	0	0	CATATTGTGTGATCACCTGTCAGCTAAGGACTCAAGACCATACCCATACTCTTCTGCTGTACTGGTTTTACCAGCACTGAGGCTTAAATAGCTAGTAATA	>CDDBBEEEECA@>?FDGHGHCE@CGCCECIGGGFF<IGGE8IJJIHFIJJIIJJJJGIHEJJJJJJIHGJIIIJHJIIIJJJJJJJHHHHHFFFFFCCC	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:100";
		String readStr = "chr11_118380033_118380233_2__numedges:341_totaledgecounts:133107_medianedgecount:417_minedgecount:4_terminatedatrepeat:false	0	chr11	118379934	239	235M3I164M	*	0	0	GGTTAGAATCAGAGAATATCAATGCTAAAAGGATTATGAGAATCACCCACTTTACCTACTTATTTTCTACATTTAAAAAAAAAAATCTAAGCTCCAAAGAAGTTAAGTGATTTGGCCCACATTGGACTGAAACTTGGCGCACCTGTCTCTCGGTGCAGTGTTCTTCCAGTACATATTGTGTGATCACCTGTCAGCTAAGGACTCAAGACCATACCCATACTCTTCTGCTGTACTGGCTGTTTTACCAGCACTGAGGCTTAAATAGCTAGTAATAACCTGACTTCACTTTTTAGTTGTTACTAAAGAAAACTAAGAACCATTTTTATTAGATAGTCAGATTTTGGTTACAATACCAGATACATCTCCATGGCATTTTCCATCAGTTCTAATGAATTTGATTAG	*	XE:i:9	XF:i:0	XN:i:0	AS:i:388	XS:i:0";
		
		SamStringReader rdr = new SamStringReader();
		
		SAMRecord read = rdr.getRead(readStr);
		
		System.out.println("read: " + read);
		
	}
	*/
}
