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
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Converts a SAM or BAM file in Genome coordinates to Transcriptome coordinates.
 * 
 * @author lmose
 */
public class GenomeToTranscriptomeConverter {
    
    private IsoformIndex isoformIndex;
    private boolean isPositiveStrandReportingOnly = true;
    private boolean shouldOutputXgTags = false;
    private IsoformOrderLoader isoformOrderLoader;
    
    private boolean isSingleEnd;
    private ReverseComplementor reverseComplementor = new ReverseComplementor();
        
    private int totalPairsOutput = 0;
    
    public GenomeToTranscriptomeConverter(IsoformIndex isoformIndex, IsoformOrderLoader isoformOrderLoader,
    		boolean isSingleEnd) {
        this.isoformIndex = isoformIndex;
        this.isoformOrderLoader = isoformOrderLoader;
        this.isSingleEnd = isSingleEnd;
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
        List<Isoform> isoforms = new ArrayList<Isoform>(isoformIndex.getAllIsoforms());
        
        // Sort the isoforms if an isoformOrderLoader was specified
        if (isoformOrderLoader != null) {
            Collections.sort(isoforms, new Isoform.IsoformOrderComparator(isoformOrderLoader));
            isoformOrderLoader.clearCache();
        }
        return isoforms;
    }
    
    private List<Isoform> getPotentialIsoforms(SAMRecord read1, SAMRecord read2) {
    	List<Isoform> isoforms = null;
    	
    	if (read1.getReferenceName().equals(read2.getReferenceName())) {
    		isoforms = getPotentialIsoforms(read1.getReferenceName(), read1.getAlignmentStart(), read2.getAlignmentEnd());
    	} else {
    		isoforms = new ArrayList<Isoform>();
    	}
    	
    	return isoforms;
    }
    
    /** 
     * Returns a list all isoforms that approximately match the specified
     * genomic coordinates. 
     */
    public List<Isoform> getPotentialIsoforms(String chromosome, int genomicStartPos, int genomicEndPos) {
        List<Isoform> potentialIsoforms = new ArrayList<Isoform>();
        Collection<String> potentials = isoformIndex.getPotentialIsoforms(chromosome, genomicStartPos);
        
        for (String isoformId : potentials) {
            Isoform isoform = isoformIndex.getIsoform(isoformId);
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
        
    private SAMRecord buildTranscriptRead(SAMRecord read, Isoform isoform, 
    		SAMFileHeader header, Cigar positiveCigar, Cigar negativeCigar, int sequenceLength, List<Coordinate> readCoordinates) {
    	
        SAMRecord transcriptRead = cloneRead(read);
        
        transcriptRead.setHeader(header);
        
        // Set ref name to the isoform
        transcriptRead.setReferenceName(isoform.getIsoformId());
                                    
        Cigar cigar = null;
        int readAlignmentStart;
        
        if ((isPositiveStrandReportingOnly) || (isoform.isPositiveStrand())) {
            cigar = positiveCigar;
            readAlignmentStart = readCoordinates.get(0).getStart();
            
        } else {
            // Negative strand
            // invert sequenceLength
            sequenceLength *= -1;
            cigar = negativeCigar;

            // Reference alignments from the end of the isoform
            readAlignmentStart = isoform.getLength() - readCoordinates.get(readCoordinates.size()-1).getStop() + 1;
            
            transcriptRead.setReadBases(reverseComplementor.reverseComplement(transcriptRead.getReadBases()));
            transcriptRead.setBaseQualities(reverseComplementor.reverse(transcriptRead.getBaseQualities()));
            
            transcriptRead.setReadNegativeStrandFlag(!read.getReadNegativeStrandFlag());
            
            // Check for single end?
            transcriptRead.setMateNegativeStrandFlag(!read.getMateNegativeStrandFlag());
        }
        
        transcriptRead.setAlignmentStart(readAlignmentStart);
        
        transcriptRead.setInferredInsertSize(sequenceLength);
        
        transcriptRead.setCigar(cigar);
                
        // For troubleshooting purposes
        if (shouldOutputXgTags) {
            transcriptRead.setAttribute("XG", read.getReferenceName() + "." + read.getAlignmentStart());
        }
        
        return transcriptRead;
    }
    
    // Single end, does not include dupe or cluster counting.
    private void convertAndOutput(SAMRecord read, SAMFileWriter outputSam, SAMFileHeader header) {
                
        List<Isoform> potentialIsoformMatches = getPotentialIsoforms(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
        
        Cigar positiveCigar = null;
        Cigar negativeCigar = null;
        
        // Convert to transcript Cigar
        if (!potentialIsoformMatches.isEmpty()) {
            List<CigarElement> elements = stripIntrons(read.getCigar());
            
            positiveCigar = getPositiveStrandCigar(elements);
            // Mutates elements!
            negativeCigar = getNegativeStrandCigar(elements);
        }
        
        for (Isoform isoform : potentialIsoformMatches) {
            List<Coordinate> readCoordinates = isoform.match(read);
            
            if (!readCoordinates.isEmpty()) {
                
                // Calc the length of the isoform insert (coordinates are inclusive so add 1 to get correct len)
                int sequenceLength = 
                    readCoordinates.get(readCoordinates.size()-1).getStop() -
                    readCoordinates.get(0).getStart() + 1;
                
                SAMRecord transcriptRead = buildTranscriptRead(read, isoform, header, positiveCigar, negativeCigar, sequenceLength, readCoordinates);
                // Mate info is unspecified
                transcriptRead.setMateAlignmentStart(0);
                transcriptRead.setMateReferenceName("*");
                transcriptRead.setMateUnmappedFlag(true);
                
                outputSam.addAlignment(transcriptRead);
            }
        }
    }

    // Paired end
    private void convertAndOutput(SAMRecord read1, SAMRecord read2, SAMFileWriter outputSam, SAMFileHeader header) {

        List<Isoform> potentialIsoformMatches = getPotentialIsoforms(read1, read2);
        
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
                        
            if (!read1Coordinates.isEmpty() && !read2Coordinates.isEmpty()) {
                
                // Calc the length of the isoform insert (coordinates are inclusive so add 1 to get correct len)
                int sequenceLength = 
                    read2Coordinates.get(read2Coordinates.size()-1).getStop() -
                    read1Coordinates.get(0).getStart() + 1;
                
                SAMRecord transcriptRead1 = buildTranscriptRead(read1, isoform, header, positiveCigar1, negativeCigar1, sequenceLength, read1Coordinates);
                SAMRecord transcriptRead2 = buildTranscriptRead(read2, isoform, header, positiveCigar2, negativeCigar2, -sequenceLength, read2Coordinates);
                
                transcriptRead1.setMateAlignmentStart(transcriptRead2.getAlignmentStart());
	            transcriptRead2.setMateAlignmentStart(transcriptRead1.getAlignmentStart());
	            transcriptRead1.setMateReferenceName(isoform.getIsoformId());
	            transcriptRead2.setMateReferenceName(isoform.getIsoformId());

                outputSam.addAlignment(transcriptRead1);
                outputSam.addAlignment(transcriptRead2);
                
                totalPairsOutput++;
            }
        }
    }
    
    private void convertFileForSingleEnd(String inputFileName, SAMFileWriter outputSam, SAMFileHeader header) {
    	System.out.println("Processing single end reads");
    	
        File inputFile = new File(inputFileName);
        
        SAMFileReader reader = new SAMFileReader(inputFile);
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        for (SAMRecord read : reader) {
        	convertAndOutput(read, outputSam, header);
        }
        
        reader.close();
    }
    
    private void convertFileForPairedEnd(String inputFileName, SAMFileWriter outputSam, SAMFileHeader header) {
    	System.out.println("Processing paired end reads");
        SamReadPairReader reader = new SamReadPairReader(inputFileName);
        for (ReadPair readPair : reader) {
            convertAndOutput(readPair.getRead1(), readPair.getRead2(), outputSam, header);
        }
        reader.close();
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
  
        if (isSingleEnd) {
        	convertFileForSingleEnd(inputFileName, outputSam, header);
        } else {
        	convertFileForPairedEnd(inputFileName, outputSam, header);
        }
        
        outputSam.close();
        
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
//        System.out.println("Starting.");
//        long s = System.currentTimeMillis();
//        
//        System.out.println("Parsing arguments");
//        Args argz = new Args(args);
//        
//        System.out.println("Determining isoform header order");
//        IsoformOrderLoader isoformOrderLoader = new IsoformOrderLoader();
//        isoformOrderLoader.loadOrdering(argz.getOrderFastaFile());
//        
//        System.out.println("Building read index");
//        BedReader bedReader = new BedReader();
//        bedReader.buildReadToIsoformIndex(argz.getBedFile(), argz.getReadOffset());
//        
//        System.out.println("Converting");
//
//        String dupeFile = argz.getDupeFile();
//        GenomeToTranscriptomeConverter converter = new GenomeToTranscriptomeConverter(bedReader, isoformOrderLoader, dupeFile);
//        converter.convertFile(argz.getInputAlignmentFile(), argz.getOutputAlignmentFile());
//                
//        long e = System.currentTimeMillis();
//        
//        System.out.println("Elapsed: " + (e-s)/1000);
    }
}
