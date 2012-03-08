package edu.unc.bioinf.ubu.sam;

import java.util.Date;

/**
 * Main class for converting genome to transcriptome coordinates
 *  
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class GenomeToTranscriptome {

    public static void main(String[] args) throws Exception {
        run(args);

    	/*
        String argz = 
            "--bed /home/lisle/gaf/ref/gaf.bed " +
            "--offset 25 " +
            "--order /home/lisle/data/coord_convert/hg19_M_ref.transcripts.fa " +
            "--xgtags " +
            "--reverse " +
            "--in /home/lisle/data/xlate/self_mate/sorted.bam " +
//            "--single " +
            "--out /home/lisle/data/xlate/self_mate/out.bam "
//            "--in /home/lisle/gtof/conv/really_small_sorted_by_ref_and_read.bam " +
//            "--out /home/lisle/gtof/conv/transcriptome_lazy.bam "
//            "--in /home/lisle/data/coord_convert/sorted_tiny.sam " +
//            "--out /home/lisle/data/coord_convert/tiny_isoform2.sam "
//          "--in /home/lisle/data/coord_convert/testing2/reverse.sam " +
//          "--out /home/lisle/data/coord_convert/testing2/converted_reverse2.bam "
//            "--in /home/lisle/data/coord_convert/testing4/gen1.sam " +
//            "--out /home/lisle/data/coord_convert/testing4/converted_gen1.sam"
//          "--in /home/lisle/data/coord_convert/testing3/sorted_genomic.sam " +
//          "--out /home/lisle/data/coord_convert/testing3/converted_sorted_genomic.bam "
            ;
        
        run(argz.split(" "));
        */

        System.out.println("mem: " + Runtime.getRuntime().totalMemory());
    }
    
    public static void run(String[] args) throws Exception {
        System.out.println("Starting " + new Date());
        long s = System.currentTimeMillis();
        
        // Parse options
        System.out.println("Parsing arguments");
        GenomeToTranscriptomeOptions options = new GenomeToTranscriptomeOptions();
        options.parseOptions(args);
        
        if (options.isValid()) {
            
            System.out.println(getAllArgsString(options));
            
            // Load isoform ordering if fasta file is specified
            IsoformOrderLoader isoformOrderLoader = null;
            
            if (options.hasOrderingFastaFile()) {
                System.out.println("Determining isoform header order");
                isoformOrderLoader = new IsoformOrderLoader();
                isoformOrderLoader.loadOrdering(options.getOrderingFastaFile());
            } else {
                System.out.println("Skipping isoform header order determination");
            }
            
            // Build read index from bed file
            System.out.println("Building read index");
            IsoformIndex isoformIndex = new IsoformIndex();
            isoformIndex.buildReadToIsoformIndex(options.getBedFile(), options.getReadOffset());
            
            // Instantiate converter and run
            System.out.println("Converting");
            
            GenomeToTranscriptomeConverter converter = new GenomeToTranscriptomeConverter(isoformIndex, isoformOrderLoader, options.isSingleEnd());
            converter.setPositiveStrandReportingOnly(options.isPositiveStrandReportingOnly());
            converter.setShouldOutputXgTags(options.shouldOutputXgTags());
            
            converter.convertFile(options.getInputAlignmentFile(), options.getOutputAlignmentFile());
                    
            long e = System.currentTimeMillis();
            
            System.out.println("Finished: " + new Date());
            System.out.println("Elapsed seconds: " + (e-s)/1000);
        }
    }
    
    private static String getAllArgsString(GenomeToTranscriptomeOptions options) {
        StringBuffer buf = new StringBuffer();
        
        buf.append("bedFile: " + options.getBedFile() + "\n");            
        buf.append("readOffset: " + options.getReadOffset() + "\n");
        buf.append("inputAlignmentFile: " + options.getInputAlignmentFile() + "\n");
        buf.append("outputAlignmentFile: " + options.getOutputAlignmentFile() + "\n");
        
        if (options.hasOrderingFastaFile()) {
            buf.append("orderFastaFile: " + options.getOrderingFastaFile() + "\n");
        } else {
            buf.append("orderFastaFile: not specified\n");
        }
                
        buf.append("positiveStrandReportingOnly: " + options.isPositiveStrandReportingOnly() + "\n");
        
        buf.append("xgtags: " + options.shouldOutputXgTags() + "\n");
        
        buf.append("single: " + options.isSingleEnd() + "\n");
        
        return buf.toString();
    }
}
