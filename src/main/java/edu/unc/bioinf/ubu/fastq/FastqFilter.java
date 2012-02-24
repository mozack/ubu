package edu.unc.bioinf.ubu.fastq;


/**
 * Filtering class for discarding unmatched reads in paired end fastq files.
 * 
 * @author lmose
 */
@Deprecated
public class FastqFilter {
    
    private int maxCachedLines = 100;
    
    private String filename1;
    private String filename2;
    private String matchFilename1;
    private String matchFilename2;
    
    private FastqInputFile file1;
    private FastqInputFile file2;
    private FastqOutputFile matchFile1;
    private FastqOutputFile matchFile2;
        
    public FastqFilter(String filename1, String filename2, String matchFilename1, String matchFilename2,
            int maxCachedLines) {
        
        this.filename1      = filename1;
        this.filename2      = filename2;
        this.matchFilename1 = matchFilename1;
        this.matchFilename2 = matchFilename2;
        this.maxCachedLines = maxCachedLines;
        
        init();
    }
    
    public String getInputFile1() {
        return filename1;
    }
    
    public String getInputFile2() {
        return filename2;
    }

    public String getMatchFile1() {
        return matchFilename1;
    }

    public String getMatchFile2() {
        return matchFilename2;
    }
    
    /**
     *  Alternate constructor to allow for dependency injection.
     *  We should look into using Spring to accomplish this. 
     */
    public void init(FastqInputFile file1, FastqInputFile file2,
            FastqOutputFile matchFile1, FastqOutputFile matchFile2) {
        
        this.file1 = file1;
        this.file2 = file2;
        this.matchFile1 = matchFile1;
        this.matchFile2 = matchFile2;
    }
    
    public void init() {
        init(new FastqInputFile(), new FastqInputFile(), new FastqOutputFile(), new FastqOutputFile());
    }
    
    public void filter() throws Exception {
        long start = System.currentTimeMillis();
        
        file1.init(filename1, maxCachedLines);
        file2.init(filename2, maxCachedLines);
        matchFile1.init(matchFilename1);
        matchFile2.init(matchFilename2);
        
        int idx1 = 1;
        int idx2 = 1;
        
        boolean isDone = false;
        
        while (!isDone) {
            FastqRecord rec1 = file1.getRecord(idx1);
            FastqRecord rec2 = file2.getRecord(idx2);
            
            if ((rec1 == null) && (rec2 == null)) {
                isDone = true;
            } else if (rec1 == null) {
                // Discard the rest of file2
                isDone = true;
            } else if (rec2 == null) {
                // Discard the rest of file1
                isDone = true;
            } else {
                if (rec1.hasSameBaseId(rec2)) {
                    idx1++;
                    idx2++;
                    matchFile1.write(rec1);
                    matchFile2.write(rec2);
                } else {
                    // Search for a match for rec1 in file2
                    boolean isFound = false;
                    int cnt = 1;
                    while ((!isFound) && (cnt < maxCachedLines) && (rec2 != null)) {
                        rec2 = file2.getRecord(idx2+cnt);
                        if (rec1.hasSameBaseId(rec2)) {
                            idx2 = idx2 + cnt + 1;
                            idx1++;
                            isFound = true;
                            matchFile1.write(rec1);
                            matchFile2.write(rec2);
                        } else {
                            cnt++;
                        }
                    }
                    
                    // If no match has been found for rec1, Search for a match for rec2 in file1
                    if (!isFound) {
                        rec2 = file2.getRecord(idx2);
                        cnt = 1;
                        while ((!isFound) && (cnt < maxCachedLines) && (rec1 != null)) {
                            rec1 = file1.getRecord(idx1+cnt);
                            if (rec2.hasSameBaseId(rec1)) {
                                idx1 = idx1 + cnt + 1;
                                idx2++;
                                isFound = true;
                                matchFile1.write(rec1);
                                matchFile2.write(rec2);
                            } else {
                                cnt++;
                            }
                        }
                    }
                    
                    if (!isFound) {
                        // The current record in both files does not match.  Move to the next 2 records.
                        idx1++;
                        idx2++;
                    }
                }
            }
            
            if ( ((idx1 % 1000000) == 0) || ((idx2 % 1000000) == 0) ) {
                System.out.println("Idx: " + idx1 + " - " + idx2);
            }
        }
        
        matchFile1.close();
        matchFile2.close();
        
        long end = System.currentTimeMillis();
        
        System.out.println("FastqFilter elapsed time(seconds): " + (end-start)/1000);
    }

    public static void main(String[] args) throws Exception {
        System.out.println("THIS CLASS IS DEPRECATED!  DO NOT USE IT!");
        
        int maxCachedLines = Integer.parseInt(args[4]);
        
        // Args: input1, input2, output1, output2, linesToCache
        FastqFilter c = new FastqFilter(args[0], args[1], args[2], args[3], maxCachedLines); 
        
        c.filter();
        
        System.out.println("Done");
        System.out.println("THIS CLASS IS DEPRECATED!  DO NOT USE IT!");
    }
}

