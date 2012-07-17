package edu.unc.bioinf.ubu.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Test utility that may be used to load text file content
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FileLoader {

    public String loadFileContent(String filename) throws IOException {
        
        StringBuffer content = new StringBuffer();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        
        String line = null;
        
        while ((line = reader.readLine()) != null) {
            content.append(line);
            content.append('\n');
        }
        
        return content.toString();
    }
}
