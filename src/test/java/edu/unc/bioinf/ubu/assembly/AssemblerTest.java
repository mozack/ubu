package edu.unc.bioinf.ubu.assembly;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.testng.Assert;

import org.testng.annotations.Test;

import edu.unc.bioinf.ubu.assembly.Assembler;

public class AssemblerTest {
	
	private static final String INPUT_FILE = "src/test/java/edu/unc/bioinf/ubu/assembly/input1.fastq";
	private static final String OUTPUT_FILE = "src/test/java/edu/unc/bioinf/ubu/assembly/output1.fasta";

	@Test( groups = "unit" )
	public void testBasicAssembly() throws IOException {
//    	File output = new File(OUTPUT_FILE);
//    	output.delete();
//		
//		Assembler ayc = new Assembler();
//		ayc.setKmerSize(3);
//		ayc.setMinNodeFrequncy(2);
//		ayc.setMinEdgeFrequency(1);
//		ayc.setMinContigLength(4);
//		ayc.setMinEdgeRatio(.1);
//		ayc.assemble(INPUT_FILE, OUTPUT_FILE);
//		
//		BufferedReader reader = new BufferedReader(new FileReader(OUTPUT_FILE));
//		String line1 = reader.readLine();
//		String line2 = null;
//		if (line1 != null) {
//			line2 = reader.readLine();
//		}
//		
//		assertTrue(line1.startsWith(">contig0"));
//		assertEquals("GCTCCAG", line2);
	}
}
