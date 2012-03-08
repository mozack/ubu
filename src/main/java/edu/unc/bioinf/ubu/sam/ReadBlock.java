package edu.unc.bioinf.ubu.sam;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Represents a Block of a read.  i.e. A Cigar of 15M5I30M would be represented
 * by 3 ReadBlocks.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class ReadBlock {
    private int readStart;
    private int referenceStart;
    private int length;
    private CigarOperator type;
    
    ReadBlock(int readStart, int referenceStart, int length, CigarOperator type) {
        this.readStart = readStart;
        this.referenceStart = referenceStart;
        this.length = length;
        this.type = type;
    }

    public int getReadStart() {
        return readStart;
    }

    public int getReferenceStart() {
        return referenceStart;
    }
    
    public int getReferenceStop() {
        return referenceStart + length - 1;
    }

    public int getLength() {
        return length;
    }

    public CigarOperator getType() {
        return type;
    }
    
    //TODO - Move elsewhere and make non-static
    public static List<ReadBlock> getReadBlocks(SAMRecord read) {    
        final Cigar cigar = read.getCigar();
        if (cigar == null) return Collections.emptyList();

        final List<ReadBlock> readBlocks = new ArrayList<ReadBlock>();
        int readBase = 1;
        int refBase  = read.getAlignmentStart();

        for (final CigarElement e : cigar.getCigarElements()) {
            
            readBlocks.add(new ReadBlock(readBase, refBase, e.getLength(), e.getOperator()));
            
            switch (e.getOperator()) {
//                case H : break; // ignore hard clips
//                case P : break; // ignore pads
                case S : readBase += e.getLength();
//                System.out.println(read.getReadName());
                break; // soft clip read bases
                case N : refBase += e.getLength(); break;  // reference skip
                case D : refBase += e.getLength(); break;
                case I : readBase += e.getLength(); break;
                case M :
//                case EQ :
//                case X :
                    final int length = e.getLength();
                    readBase += length;
                    refBase  += length;
                    break;
                default : throw new IllegalStateException(
                        "Case statement didn't deal with cigar op: " + e.getOperator() +
                        " for read: [" + read.getReadName() + "]");
            }
        }
        
        return Collections.unmodifiableList(readBlocks);
    }
}
