package edu.unc.bioinf.ubu.sam;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Builds a matrix of expected counts from rsem files.
 *  
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class RsemExpectedCountsAggregator {
    
    private List<RsemFileInfo> fileInfo = new ArrayList<RsemFileInfo>();
    private BufferedWriter writer;
    
    private List<BufferedReader> readers = new ArrayList<BufferedReader>(); 
    
    private int lineNum = 1;
    
    private String countFile;

    public void aggregate(String countFile, String outputFile, List<String> directories) throws IOException {
        this.countFile = countFile;
        writer = new BufferedWriter(new FileWriter(outputFile, false));
        findRsemFiles(directories);
        buildMatrix();
        writer.close();
        closeReaders();
    }
    
    private void closeReaders() throws IOException {
        for (BufferedReader reader : readers) {
            reader.close();
        }
    }
    
    private void buildMatrix() throws IOException {
        initReaders();
        writeHeader();
        processFiles();
    }
    
    private List<String> getNextResultList() throws IOException {
        List<String> results = new ArrayList<String>();
        boolean isNullInResults = false;
        boolean isNonNullInResults = false;
        
        for (BufferedReader reader : readers) {
            String line = reader.readLine();
            if (line == null) {
                isNullInResults = true;
            }
            
            if (line != null) {
                isNonNullInResults = true;
            }
            results.add(line);
        }
        
        if (isNullInResults && isNonNullInResults) {
            throw new IllegalArgumentException("Invalid EOF reached at line: " + lineNum);
        }
        
        lineNum++;
        
        if (isNullInResults) {
            results = null;
        }
        
        return results;
    }
    
    private void processFiles() throws IOException {
        List<String> lines = getNextResultList();

        while (lines != null) {
            if (containsNonZeroCount(lines)) {
                outputCount(lines);
            }
            
            lines = getNextResultList();
        }
    }
    
    private boolean containsNonZeroCount(List<String> lines) {
        double totalCount = 0.0;
        
        for (String line : lines) {
            String[] fields = line.split("\t");
            String countStr = fields[1];
            
            totalCount += Double.valueOf(countStr);
        }
        
        boolean isNonZero = false;
        
        if (Math.abs(totalCount) > .0001) {
            isNonZero = true;
        }
        
        return isNonZero;
    }
    
    private void outputCount(List<String> lines) throws IOException {
        StringBuffer geneLine = new StringBuffer();
        String firstLine = lines.get(0);
        String[] fields = firstLine.split("\t");
        String geneId = fields[0];
        geneLine.append(geneId);
        
        for (String line : lines) {
            String[] lineFields = line.split("\t");
            String currGeneId = lineFields[0];
            if (!(currGeneId.equals(geneId))) {
                throw new IllegalArgumentException ("Mismatch genes: " + geneId + " - " + currGeneId);
            }
            
            String expectedCount = lineFields[1];
            
            geneLine.append('\t');
            geneLine.append(expectedCount);
        }
        
        writer.write(geneLine.toString());
        writer.write('\n');
    }
    
    private void initReaders() throws IOException {
        for (RsemFileInfo info : fileInfo) {
            readers.add(new BufferedReader(new FileReader(info.getGeneFile())));
        }
    }
    
    private void writeHeader() throws IOException {
        StringBuffer str = new StringBuffer();
        for (RsemFileInfo info : fileInfo) {
            // Skip first column
            str.append('\t');
            str.append(info.getSampleId());
        }
        
        writer.write(str.toString());
        writer.write('\n');
    }
    
    private void findRsemFiles(List<String> directories) {
        for (String directory : directories) {
            System.out.println("Processing directory: " + directory);
            File file = new File(directory);
            
            File[] filesInDir = file.listFiles();
            for (File fileInDir : filesInDir) {
                if (fileInDir.isDirectory()) {
                    System.out.println("Adding: " + fileInDir.getName());
                    fileInfo.add(new RsemFileInfo(
                            fileInDir.getName(),
                            fileInDir.getAbsolutePath() + "/" + countFile,
                            fileInDir.getAbsolutePath() + "/" + "foo"));
                }
            }
        }
    }
    
    static class RsemFileInfo {
        private String sampleId;
        private String geneFile;
        private String isoformFile;
        
        RsemFileInfo(String sampleId, String geneFile, String isoformFile) {
            this.sampleId = sampleId;
            this.geneFile = geneFile;
            this.isoformFile = isoformFile;
        }
        
        public String getSampleId() {
            return sampleId;
        }
        
        public String getGeneFile() {
            return geneFile;
        }
        
        public String getIsoformFile() {
            return isoformFile;
        }
    }
    
    public static void main(String[] args) throws IOException {
        
        if (args.length < 3) {
            System.out.println("RsemExpectedCountsAggregator <results_file_name> <output_file> <dir1> <dir2> ... <dirn>");
            System.exit(0);
        }

        
//        args = new String[] {
//                "rsem.genes.results",
//                "all_gene_counts2.tsv",
//                "/home/lisle/data/rsem"
//        };

        
        String countFile = args[0];
        String outputFile = args[1];
        
        List<String> directories = new ArrayList<String>();
        
        for (int i=2; i<args.length; i++) {
            directories.add(args[i]);
        }
        
        RsemExpectedCountsAggregator aggregator = new RsemExpectedCountsAggregator();
        aggregator.aggregate(countFile, outputFile, directories);
    }
}

