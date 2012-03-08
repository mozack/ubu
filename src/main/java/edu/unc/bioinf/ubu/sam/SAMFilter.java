package edu.unc.bioinf.ubu.sam;

import java.io.File;
import java.io.IOException;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Filters paired reads from a SAM or BAM file.  Indels and/or clusters
 * greater than a specified insert length may be filtered.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SAMFilter {
    
    private static final String INPUT_FILE = "input";
    private static final String OUTPUT_FILE = "output";
    private static final String STRIP_INDELS = "strip-indels";
    private static final String MAX_INSERT_LEN = "max-insert";
    
    private boolean shouldStripIndels = false;
    private int     maxInsertLen = -1;

    public void filter(String input, String output) {
        File outputFile = new File(output);
        SamReadPairReader reader = new SamReadPairReader(input);
        
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getHeader(),
                false, outputFile);
        
        for (ReadPair pair : reader) {
            SAMRecord read1 = pair.getRead1();
            SAMRecord read2 = pair.getRead2();
            
            boolean isPairIncluded = !shouldStripIndels || !hasIndels(read1, read2);
            
            if (isPairIncluded) {
                isPairIncluded = !isMaxInsertLenExceeded(read1, read2);
            }
            
            // Output only the pairs that have passed our tests
            if (isPairIncluded) {
                writer.addAlignment(read1);
                writer.addAlignment(read2);
            }
        }
        
        writer.close();
        reader.close();
    }
    
    private boolean isMaxInsertLenExceeded(SAMRecord read1, SAMRecord read2) {
        boolean isMaxExceeded = false;
        
        if (isMaxInsertLenSpecified()) {
            isMaxExceeded = 
                (Math.abs(read1.getInferredInsertSize()) > maxInsertLen) &&
                (Math.abs(read2.getInferredInsertSize()) > maxInsertLen);
        }
        
        return isMaxExceeded;
    }
    
    private void setMaxInsertLen(int len) {
        maxInsertLen = len;
    }
    
    private boolean isMaxInsertLenSpecified() {
        return maxInsertLen > 0;
    }
    
    private boolean hasIndels(SAMRecord read1, SAMRecord read2) {
        return ( (hasIndel(read1)) || (hasIndel(read2)) );
    }
    
    public void setShouldStripIndels(boolean shouldStripIndels) {
        this.shouldStripIndels = true;
    }
    
    private boolean hasIndel(SAMRecord read) {
        for (CigarElement element : read.getCigar().getCigarElements()) {
            if ((element.getOperator() == CigarOperator.D) ||
                (element.getOperator() == CigarOperator.I)) {
                return true;
            }
        }
        
        return false;
    }
    
    private static boolean isOptionSetValid(OptionSet options, OptionParser parser) throws IOException {
        if ((!options.hasArgument(INPUT_FILE)) ||
            (!options.hasArgument(OUTPUT_FILE))) {
            
            parser.printHelpOn(System.err);
            
            return false;
        }
        
        return true;
    }
        
    public static void main(String[] args) throws IOException {
        long s = System.currentTimeMillis();
        
//        args = "--input /home/lisle/data/coord_convert/converted_small_sg.bam --output /home/lisle/data/coord_convert/stripped_converted_small_sg.bam --strip-indels --max-insert 4000".split(" ");

        OptionParser parser = new OptionParser();
        parser.accepts(INPUT_FILE, "Required input sam or bam file").withRequiredArg().ofType(String.class);
        parser.accepts(OUTPUT_FILE, "Required output sam or bam file").withRequiredArg().ofType(String.class);
        parser.accepts(STRIP_INDELS, "If specified, discard read pairs containing indels from output (default off)");
        parser.accepts(MAX_INSERT_LEN, "If specified, discard clusters greater than specified insert length").withRequiredArg().ofType(Integer.class);
        
        OptionSet options = parser.parse(args);
        
        if (isOptionSetValid(options, parser)) {
            
            SAMFilter filter = new SAMFilter();
            if (options.has(STRIP_INDELS)) {
                filter.setShouldStripIndels(true);
            }
            
            if (options.hasArgument(MAX_INSERT_LEN)) {
                filter.setMaxInsertLen((Integer) options.valueOf(MAX_INSERT_LEN));
            }
            
            filter.filter((String) options.valueOf(INPUT_FILE), (String) options.valueOf(OUTPUT_FILE));
            
    //        new IndelStripper().strip("/home/lisle/data/coord_convert/million_isoform.sam",
    //                "/home/lisle/data/coord_convert/stripped_million_isoform.sam");
            
            long e = System.currentTimeMillis();
            
            System.out.println("Elapsed: " + (e-s)/1000);
        }
    }
}
