package edu.unc.bioinf.ubu.assembly;

public class SequenceUtil {

	public static boolean isMatch(String seq1, String seq2, int allowedMismatches) {
		
		if (allowedMismatches == 0) {
			return seq1.equals(seq2);
		}
		
		int mismatches = 0;
		
		if (seq1.length() != seq2.length()) {
			return false;
		}
		
		for (int i=0; i<seq1.length(); i++) {
			if (seq1.charAt(i) != seq2.charAt(i)) {
				mismatches += 1;
			}
			
			if (mismatches > allowedMismatches) {
				return false;
			}
		}
		
		return true;
	}
}
