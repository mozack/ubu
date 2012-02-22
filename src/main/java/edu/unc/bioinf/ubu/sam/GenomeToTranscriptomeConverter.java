package edu.unc.bioinf.ubu.sam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Converts a SAM or BAM file in Genome coordinates to Transcriptome coordinates.
 * 
 * @author lmose
 */
public class GenomeToTranscriptomeConverter {
    
    private BedReader bedReader;
    private boolean isPositiveStrandReportingOnly = true;
    private boolean shouldOutputXgTags = false;
    private IsoformOrderLoader isoformOrderLoader;
    private String currentCluster = "";
    private int currentClusterCount = 0;
    // key = dupeCount, value = # clusters w/ this dupe count
    private Map<Integer, Integer> dupeFrequencyMap = new HashMap<Integer, Integer>();
    
    private String dupeFile;
    private ReverseComplementor reverseComplementor = new ReverseComplementor();
    
    private boolean is1ReadMappedForCluster = false;
    private boolean isPairMappedForCluster = false;
    
    private int totalClusters = 0;
    private int oneReadClusters = 0;
    private int pairedClusters = 0;
    
    private int totalPairsOutput = 0;
    
    public GenomeToTranscriptomeConverter(BedReader bedReader, IsoformOrderLoader isoformOrderLoader, String dupeFile) {
        this.bedReader = bedReader;
        this.isoformOrderLoader = isoformOrderLoader;
        this.dupeFile = dupeFile;
    }
   
    private SAMFileHeader buildHeader() {
        SAMFileHeader header = new SAMFileHeader();
        
        for (Isoform isoform : getSortedIsoforms()) {
            SAMSequenceRecord sq = new SAMSequenceRecord(isoform.getIsoformId(), isoform.getLength());
            header.addSequence(sq);            
        }
        
        return header;
    }
    
    /**
     * Returns a list of all Isoforms available sorted by id.
     */
    private List<Isoform> getSortedIsoforms() {
        List<Isoform> isoforms = new ArrayList<Isoform>(bedReader.getAllIsoforms());
        
        // Sort the isoforms if an isoformOrderLoader was specified
        if (isoformOrderLoader != null) {
            Collections.sort(isoforms, new Isoform.IsoformOrderComparator(isoformOrderLoader));
            isoformOrderLoader.clearCache();
        }
        return isoforms;
    }
    
    private boolean shouldTrackDupes() {
        return dupeFile != null;
    }
    
    /** 
     * Returns a list all isoforms that approximately match the specified
     * genomic coordinates. 
     */
    public List<Isoform> getPotentialIsoforms(String chromosome, int genomicStartPos, int genomicEndPos) {
        List<Isoform> potentialIsoforms = new ArrayList<Isoform>();
        Collection<String> potentials = bedReader.getPotentialIsoforms(chromosome, genomicStartPos);
        
        for (String isoformId : potentials) {
            Isoform isoform = bedReader.getIsoform(isoformId);
            if (isoform.containsWithinGenomicRange(genomicStartPos, genomicEndPos)) {
                potentialIsoforms.add(isoform);
            }
        }
        
        return potentialIsoforms;
    }
    
    private List<CigarElement> stripIntrons(Cigar cigar) {
        List<CigarElement> elements = new ArrayList<CigarElement>();
        CigarElement prevNonIntronElement = null;
        
        for (CigarElement element : cigar.getCigarElements()) {
            // Omit introns
            if (element.getOperator() != CigarOperator.N) {
                
                // Combine contiguous elements with the same operator
                if ((prevNonIntronElement != null) && (prevNonIntronElement.getOperator() == element.getOperator())) {
                    elements.remove(prevNonIntronElement);
                    element = new CigarElement(prevNonIntronElement.getLength() + element.getLength(), element.getOperator());
                    elements.add(element);
                } else {                
                    elements.add(element);
                }
                
                prevNonIntronElement = element;
            }
        }
        
        return elements;        
    }
    
    private Cigar getPositiveStrandCigar(List<CigarElement> elements) {
        Cigar cigar = new Cigar();
        
        for (CigarElement element : elements) {
            cigar.add(element);
        }
        
        return cigar;
    }
    
    private Cigar getNegativeStrandCigar(List<CigarElement> elements) {
        Collections.reverse(elements);
        Cigar cigar = new Cigar();
        
        for (CigarElement element : elements) {
            cigar.add(element);
        }
        
        return cigar;
    }
    
    private SAMRecord cloneRead(SAMRecord read) {
        try {
            return (SAMRecord) read.clone();
        } catch (CloneNotSupportedException e) {
            // Infamous "this should never happen" comment here.
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
    
    //TODO - Dupe of SamReadPairReader logic
    private String getBaseName(SAMRecord read) {
        return read.getReadName().substring(0, read.getReadName().length()-2);
    }
    
    private void updateDupeCounts() {
        updateDupes(currentClusterCount, dupeFrequencyMap);
    }
    
    private void updateDupes(int count, Map<Integer, Integer> dupeFreqMap) {
        if (count != 0) {
            Integer frequency = dupeFreqMap.get(count);
            if (frequency == null) {
                dupeFreqMap.put(count, 1);
            } else {
                dupeFreqMap.put(count, frequency + 1);
            }
        }        
    }
        
    private void writeDupeCounts(Map<Integer, Integer> dupeFreqMap, String file) throws IOException {
        
        List<Integer> dupes = new ArrayList<Integer>(dupeFreqMap.keySet());
        Collections.sort(dupes);
        
        BufferedWriter writer = new BufferedWriter(new FileWriter(file, false));
        
        for (Integer dupe : dupes) {
            writer.write(dupe + "," + dupeFreqMap.get(dupe));
            writer.write('\n');
        }
        
        writer.close();
    }
    
    private void writeDupeCounts() throws IOException {
        writeDupeCounts(dupeFrequencyMap, dupeFile);
    }
    
    private void updateClusterCounts() {
        totalClusters += 1;
        
        if (is1ReadMappedForCluster && !isPairMappedForCluster) {
            oneReadClusters += 1;
        }
        
        if (isPairMappedForCluster) {
            pairedClusters += 1;
        }
        
        is1ReadMappedForCluster = false;
        isPairMappedForCluster = false;
    }
    
    private void convertAndOutput(SAMRecord read1, SAMRecord read2, SAMFileWriter outputSam, SAMFileHeader header) {
        
        String cluster = getBaseName(read1);
        if (!cluster.equals(currentCluster)) {
            if (shouldTrackDupes()) {
                updateDupeCounts();
            }
            
            currentCluster = cluster;
            currentClusterCount = 0;
            
            updateClusterCounts();
        }
        
        List<Isoform> potentialIsoformMatches = getPotentialIsoforms(read1.getReferenceName(), read1.getAlignmentStart(), read2.getAlignmentEnd());
        
        Cigar positiveCigar1 = null;
        Cigar negativeCigar1 = null;
        Cigar positiveCigar2 = null;
        Cigar negativeCigar2 = null;
        
        // Convert to transcript Cigar
        if (!potentialIsoformMatches.isEmpty()) {
            List<CigarElement> elements1 = stripIntrons(read1.getCigar());
            List<CigarElement> elements2 = stripIntrons(read2.getCigar());
            
            positiveCigar1 = getPositiveStrandCigar(elements1);
            // Mutates elements!
            negativeCigar1 = getNegativeStrandCigar(elements1);
            positiveCigar2 = getPositiveStrandCigar(elements2);
            // Mutates elements!
            negativeCigar2 = getNegativeStrandCigar(elements2);
        }
        
        for (Isoform isoform : potentialIsoformMatches) {
            List<Coordinate> read1Coordinates = isoform.match(read1);
            List<Coordinate> read2Coordinates = isoform.match(read2);
            
            is1ReadMappedForCluster = !read1Coordinates.isEmpty() || !read2Coordinates.isEmpty();
            
            if (!read1Coordinates.isEmpty() && !read2Coordinates.isEmpty()) {
                isPairMappedForCluster = true;
                
                SAMRecord transcriptRead1 = cloneRead(read1);
                SAMRecord transcriptRead2 = cloneRead(read2);
                
                transcriptRead1.setHeader(header);
                transcriptRead2.setHeader(header);
                
                // Set ref name to the isoform
                transcriptRead1.setReferenceName(isoform.getIsoformId());
                transcriptRead2.setReferenceName(isoform.getIsoformId());
                                    
                // Calc the length of the isoform insert (coordinates are inclusive so add 1 to get correct len)
                int sequenceLength = 
                    read2Coordinates.get(read2Coordinates.size()-1).getStop() -
                    read1Coordinates.get(0).getStart() + 1;
                
                Cigar cigar1 = null;
                Cigar cigar2 = null;
                int read1AlignmentStart;
                int read2AlignmentStart;
                
                if ((isPositiveStrandReportingOnly) || (isoform.isPositiveStrand())) {
                    cigar1 = positiveCigar1;
                    cigar2 = positiveCigar2;
                    read1AlignmentStart = read1Coordinates.get(0).getStart();
                    read2AlignmentStart = read2Coordinates.get(0).getStart();
                    
                } else {
                    // Negative strand
                    // invert sequenceLength
                    sequenceLength *= -1;
                    cigar1 = negativeCigar1;
                    cigar2 = negativeCigar2;
                    // Reference alignments from the end of the isoform
                    read1AlignmentStart = isoform.getLength() - read1Coordinates.get(read1Coordinates.size()-1).getStop() + 1;
                    read2AlignmentStart = isoform.getLength() - read2Coordinates.get(read2Coordinates.size()-1).getStop() + 1;
                    
                    transcriptRead1.setReadBases(reverseComplementor.reverseComplement(transcriptRead1.getReadBases()));
                    transcriptRead2.setReadBases(reverseComplementor.reverseComplement(transcriptRead2.getReadBases()));

                    transcriptRead1.setBaseQualities(reverseComplementor.reverse(transcriptRead1.getBaseQualities()));
                    transcriptRead2.setBaseQualities(reverseComplementor.reverse(transcriptRead2.getBaseQualities()));
                    
                    transcriptRead1.setReadNegativeStrandFlag(!transcriptRead1.getReadNegativeStrandFlag());
                    transcriptRead1.setMateNegativeStrandFlag(!transcriptRead1.getMateNegativeStrandFlag());
                    transcriptRead2.setReadNegativeStrandFlag(!transcriptRead2.getReadNegativeStrandFlag());
                    transcriptRead2.setMateNegativeStrandFlag(!transcriptRead2.getMateNegativeStrandFlag());
                }
                
                transcriptRead1.setAlignmentStart(read1AlignmentStart);
                transcriptRead2.setAlignmentStart(read2AlignmentStart);
                
                transcriptRead1.setInferredInsertSize(sequenceLength);
                transcriptRead2.setInferredInsertSize(-sequenceLength);
                
                transcriptRead1.setCigar(cigar1);
                transcriptRead2.setCigar(cigar2);
                
                // Set mate info
                transcriptRead1.setMateAlignmentStart(transcriptRead2.getAlignmentStart());
                transcriptRead2.setMateAlignmentStart(transcriptRead1.getAlignmentStart());
                transcriptRead1.setMateReferenceName(isoform.getIsoformId());
                transcriptRead2.setMateReferenceName(isoform.getIsoformId());
                
                // For troubleshooting purposes
                if (shouldOutputXgTags) {
                    transcriptRead1.setAttribute("XG", read1.getReferenceName() + "." + read1.getAlignmentStart());
                    transcriptRead2.setAttribute("XG", read2.getReferenceName() + "." + read2.getAlignmentStart());
                }
                
                outputSam.addAlignment(transcriptRead1);
                outputSam.addAlignment(transcriptRead2);
                
                totalPairsOutput++;
                currentClusterCount++;
            }
        }
    }
    
    /**
     * Opens the file specified by inputFileName, processes all paired sam
     * records, converts from genome to transcriptome coordinates and writes
     * the output to outputFileName.<p>
     * Both SAM and BAM files are supported.  The file type is determined by
     * the file suffix.
     */
    public void convertFile(String inputFileName, String outputFileName) throws Exception {
        
        System.out.println("Reading file: " + inputFileName);
        System.out.println("Writing to file: " + outputFileName);
        
        File outputFile = new File(outputFileName);
        
        System.out.println("Writing header");
        SAMFileHeader header = buildHeader();
        final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(header,
                  true, outputFile);
  
        System.out.println("Processing reads");
        SamReadPairReader reader = new SamReadPairReader(inputFileName);
        for (ReadPair readPair : reader) {
            convertAndOutput(readPair.getRead1(), readPair.getRead2(), outputSam, header);
        }
        
        // Call this one more time to include the last cluster
        updateDupeCounts();
        updateClusterCounts();
        
        outputSam.close();
        reader.close();
        
        if (shouldTrackDupes()) {
            System.out.println("Outputing dupe counts to file: " + this.dupeFile);
            writeDupeCounts();
        }
        
        System.out.println("Total unique clusters input: " + totalClusters);
        System.out.println("Total unique clusters translated: " + pairedClusters);
        System.out.println("Total unique clusters with only 1 read matched: " + oneReadClusters);
        System.out.println("Total pairs output: " + totalPairsOutput);
        
        System.out.println("Done.");
    }
    
    public void setPositiveStrandReportingOnly(boolean isPositiveStrandReportingOnly) {
        this.isPositiveStrandReportingOnly = isPositiveStrandReportingOnly;
    }
    
    public void setShouldOutputXgTags(boolean shouldOutputXgTags) {
        this.shouldOutputXgTags = shouldOutputXgTags;
    }
    
    static class Args {
        
        private String orderFastaFile;
        private String bedFile;
        private String dupeFile;
        private String inputAlignmentFile;
        private String outputAlignmentFile;
        private int    readOffset;
        
        public Args(String[] args) {
            bedFile             = args[0];
            readOffset          = Integer.parseInt(args[1]);
            inputAlignmentFile  = args[2];
            outputAlignmentFile = args[3];
            orderFastaFile      = args[4];
            dupeFile            = args[5];
            
            
//            orderFastaFile = "/home/lisle/data/coord_convert/hg19_M_ref.transcripts.fa";
            
            //bedFile = "/home/lisle/data/coord_convert/chr16.bed";
//            bedFile = "/home/lisle/data/coord_convert/ucsc_known_gene_bed.txt";
//            readOffset = 300;
//            readOffset = 25;
//            dupeFile = "/home/lisle/data/coord_convert/dupes.txt";
            
//            inputAlignmentFile = "/home/lisle/data/coord_convert/sorted_tiny.sam";
//            outputAlignmentFile = "/home/lisle/data/coord_convert/tiny_isoform.sam";
            
//            inputAlignmentFile = "/home/lisle/data/coord_convert/repeating_input.sam";
//            outputAlignmentFile = "/home/lisle/data/coord_convert/converted_repeating_input.sam";
            
//            inputAlignmentFile = "/home/lisle/data/coord_convert/sorted_million.sam";
//            outputAlignmentFile = "/home/lisle/data/coord_convert/million_isoform.sam";
            
//              inputAlignmentFile = "/home/lisle/data/coord_convert/2/64.2.sam";
//              outputAlignmentFile = "/home/lisle/data/coord_convert/2/converted_64.2.sam";

//              inputAlignmentFile = "/home/lisle/data/coord_convert/2/UNC11-SN627_70:1:1101:15855:142198.sam";
//              outputAlignmentFile = "/home/lisle/data/coord_convert/2/converted_UNC11-SN627_70:1:1101:15855:142198.sam";


            
//          converter.convertFile("/home/lisle/data/coord_convert/sorted_million.sam", "/home/lisle/data/coord_convert/million_isoform.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/intron.sam", "/home/lisle/data/coord_convert/converted_intron.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/insert.sam", "/home/lisle/data/coord_convert/converted_insert.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/delete.sam", "/home/lisle/data/coord_convert/converted_delete.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/sorted_million.sam", "/home/lisle/data/coord_convert/million_isoform3.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/small_genomic_bwa.sam", "/home/lisle/data/coord_convert/converted_small_genomic_bwa.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/bwa1.sam", "/home/lisle/data/coord_convert/converted_bwa1.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/bwa2.sam", "/home/lisle/data/coord_convert/converted_bwa2.sam");
//          converter.convertFile("/home/lisle/data/coord_convert/bwa2.sam", "/home/lisle/data/coord_convert/converted_bwa2.bam");
//          converter.convertFile("/home/lisle/data/coord_convert/1_bwa_genomic_pair.sam", "/home/lisle/data/coord_convert/converted_1_bwa_genomic_pair.bam");
//          converter.convertFile("/home/lisle/data/coord_convert/gm/gm4.sam", "/home/lisle/data/coord_convert/gm/converted_gm4.sam");
            
            System.out.println(this);
        }

        public String getOrderFastaFile() {
            return orderFastaFile;
        }

        public String getBedFile() {
            return bedFile;
        }

        public String getDupeFile() {
            return dupeFile;
        }

        public String getInputAlignmentFile() {
            return inputAlignmentFile;
        }

        public String getOutputAlignmentFile() {
            return outputAlignmentFile;
        }

        public int getReadOffset() {
            return readOffset;
        }
        
        public String toString() {
            StringBuffer buf = new StringBuffer();
            
            buf.append("bedFile: " + bedFile + "\n");            
            buf.append("readOffset: " + readOffset + "\n");
            buf.append("inputAlignmentFile: " + inputAlignmentFile + "\n");
            buf.append("outputAlignmentFile: " + outputAlignmentFile + "\n");
            buf.append("orderFastaFile: " + orderFastaFile + "\n");
            buf.append("dupeFile: " + dupeFile + "\n");
            
            return buf.toString();
        }
    }
    
    public static void main(String[] args) throws Exception {
        System.out.println("Starting.");
        long s = System.currentTimeMillis();
        
        System.out.println("Parsing arguments");
        Args argz = new Args(args);
        
        System.out.println("Determining isoform header order");
        IsoformOrderLoader isoformOrderLoader = new IsoformOrderLoader();
        isoformOrderLoader.loadOrdering(argz.getOrderFastaFile());
        
        System.out.println("Building read index");
        BedReader bedReader = new BedReader();
        bedReader.buildReadToIsoformIndex(argz.getBedFile(), argz.getReadOffset());
        
        System.out.println("Converting");

        String dupeFile = argz.getDupeFile();
        GenomeToTranscriptomeConverter converter = new GenomeToTranscriptomeConverter(bedReader, isoformOrderLoader, dupeFile);
        converter.convertFile(argz.getInputAlignmentFile(), argz.getOutputAlignmentFile());
                
        long e = System.currentTimeMillis();
        
        System.out.println("Elapsed: " + (e-s)/1000);
    }
}
