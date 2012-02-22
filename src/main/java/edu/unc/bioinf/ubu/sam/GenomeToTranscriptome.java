package edu.unc.bioinf.ubu.sam;

import java.util.Date;

/**
 * Main class for converting genome to transcriptome coordinates
 *  
 * @author lmose
 */
public class GenomeToTranscriptome {

    public static void main(String[] args) throws Exception {
//        run(args);
        
        
        String argz = 
            "--bed /home/lisle/data/coord_convert/ucsc_known_gene_bed.txt " +
            "--dupes /home/lisle/data/coord_convert/dupes3.txt " +
            "--offset 300 " +
            "--order /home/lisle/data/coord_convert/hg19_M_ref.transcripts.fa " +
            "--xgtags " +
            "--reverse " +
//            "--in /home/lisle/data/coord_convert/sorted_tiny.sam " +
//            "--out /home/lisle/data/coord_convert/tiny_isoform2.sam "
//          "--in /home/lisle/data/coord_convert/testing2/reverse.sam " +
//          "--out /home/lisle/data/coord_convert/testing2/converted_reverse2.bam "
            "--in /home/lisle/data/coord_convert/testing4/gen1.sam " +
            "--out /home/lisle/data/coord_convert/testing4/converted_gen1.sam"
//          "--in /home/lisle/data/coord_convert/testing3/sorted_genomic.sam " +
//          "--out /home/lisle/data/coord_convert/testing3/converted_sorted_genomic.bam "

            ;
        
        
        run(argz.split(" "));

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
            BedReader bedReader = new BedReader();
            bedReader.buildReadToIsoformIndex(options.getBedFile(), options.getReadOffset());
            
            // Instantiate converter and run
            System.out.println("Converting");
            String dupeFile = options.hasDuplicatesFile() ? options.getDuplicatesFile() : null;
            GenomeToTranscriptomeConverter converter = new GenomeToTranscriptomeConverter(bedReader, isoformOrderLoader, dupeFile);
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
        
        if (options.hasDuplicatesFile()) {
            buf.append("dupeFile: " + options.getDuplicatesFile() + "\n");
        } else {
            buf.append("dupeFile: not specified\n");
        }
        
        buf.append("positiveStrandReportingOnly: " + options.isPositiveStrandReportingOnly() + "\n");
        
        buf.append("xgtags: " + options.shouldOutputXgTags() + "\n");
        
        return buf.toString();
    }
}
