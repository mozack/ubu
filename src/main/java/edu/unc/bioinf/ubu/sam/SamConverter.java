package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.io.IOException;

import edu.unc.bioinf.ubu.util.QualityConverter;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Converts SAM / BAM file content
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamConverter {

	private QualityConverter qualityConverter = new QualityConverter();
	
	public void convert(String inputFile, String outputFile) throws IOException {
        long start = System.currentTimeMillis();
        
        File inFile = new File(inputFile);
                
        SAMFileReader inputSam = new SAMFileReader(inFile);
        inputSam.setValidationStringency(ValidationStringency.SILENT);

        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
                true, new File(outputFile));
        
        int count = 0;
        
        for (SAMRecord read : inputSam) {
        	
        	read.setBaseQualityString(qualityConverter.phred64ToPhred33(read.getBaseQualityString()));
        	writer.addAlignment(read);
        	
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
        }
        
        writer.close();

        long stop = System.currentTimeMillis();
        
        System.out.println("Done.  Elapsed secs: " + (stop-start)/1000);
	}
	
	public static void run(String[] args) throws IOException {
		SamConverterOptions options = new SamConverterOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			SamConverter converter = new SamConverter();
			
			converter.convert(options.getInputFile(), options.getOutputFile());
		}
	}
}
