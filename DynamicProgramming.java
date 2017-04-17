import java.util.*;

public class DynamicProgramming {

	public static void main(String[] args) {
		
	}
	
	//Create a scoring matrix based on an alphabet and designated scores
	public HashMap<Character, HashMap<Character, Integer>> buildScoringMatrix(String alphabet, int diagScore, int offDiagScore, int dashScore) {
		int size = alphabet.length() + 1;
		HashMap<Character, HashMap<Character, Integer>> rowsMatrix = new HashMap<Character, HashMap<Character, Integer>>(); //2D HashMap
		
		char[] copy = new char[alphabet.length() + 1]; //Creates an array storing the alphabet and one empty space
		alphabet.getChars(0, size - 2, copy, 0);
		copy[size - 1] = '_';
		
		for (char r : copy) {
			HashMap<Character, Integer> colsMatrix = new HashMap<Character, Integer>();
			if (r == '_') {
				for (char c : copy) {
					colsMatrix.put(c, dashScore);
				}
			} else {
				for (char c : copy) {
					if (r == c) {
						colsMatrix.put(c, diagScore);
					} else if (c == '_') {
						colsMatrix.put(c, dashScore);
					} else {
						colsMatrix.put(c, offDiagScore);
					}
				}
			}
			rowsMatrix.put(r, colsMatrix);
		}
		return rowsMatrix;
	}
	
	//Help compute global alignment matrix
	public ArrayList<ArrayList<Integer>> globalHelper(String seqX, String seqY, HashMap<Character, HashMap<Character, Integer>> scoringMatrix) {
		int numRows = seqX.length();
		int numCols = seqY.length();
		ArrayList<ArrayList<Integer>> alignmentMatrix = new ArrayList<ArrayList<Integer>>();
		
		for (int i = 0; i < numRows; i++) {
			alignmentMatrix.add(new ArrayList<Integer>(alignmentMatrix.get(i).get(0) + scoringMatrix.get(seqX.charAt(i)).get('_')));	
		}
		for (int j = 0; j < numCols; j++) {
			alignmentMatrix.get(0).add(alignmentMatrix.get(0).get(j) + scoringMatrix.get('_').get(seqY.charAt(j)));
		}
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				int[] scoreSet = new int[3];
				scoreSet[0] = alignmentMatrix.get(i).get(j) + scoringMatrix.get(seqX.charAt(i)).get(seqY.charAt(j));
				scoreSet[1] = alignmentMatrix.get(i).get(j + 1) + scoringMatrix.get(seqX.charAt(i)).get('_');
				scoreSet[2] = alignmentMatrix.get(i + 1).get(j) + scoringMatrix.get(seqY.charAt(j)).get('_');
				int maxVal = -(9999);
				for (int k : scoreSet) {
					if (k > maxVal) {
						maxVal = k;
					}
				}
				alignmentMatrix.get(i + 1).add(maxVal);
			}
		}
		return alignmentMatrix;
	}
	
	//Compute local alignment matrix
	public ArrayList<ArrayList<Integer>> computeAlignmentMatrix(String seqX, String seqY, 
			HashMap<Character, HashMap<Character, Integer>> scoringMatrix, boolean globalFlag) {
		
		if (globalFlag) {
			return globalHelper(seqX, seqY, scoringMatrix);
		}
		
		int numRows = seqX.length();
		int numCols = seqY.length();
		ArrayList<ArrayList<Integer>> alignmentMatrix = new ArrayList<ArrayList<Integer>>();
		
		for (int i = 0; i < numRows; i++) {
			if (alignmentMatrix.get(i).get(0) + scoringMatrix.get(seqX.charAt(i)).get('_') > 0) {
				alignmentMatrix.add(new ArrayList<Integer>(alignmentMatrix.get(i).get(0) + scoringMatrix.get(seqX.charAt(i)).get('_')));
			} else {
				alignmentMatrix.add(new ArrayList<Integer>(0));
			}
		}
		for (int j = 0; j < numCols; j++) {
			if (alignmentMatrix.get(0).get(j) + scoringMatrix.get('_').get(seqY.charAt(j)) > 0) {
				alignmentMatrix.get(0).add(alignmentMatrix.get(0).get(j) + scoringMatrix.get('_').get(seqY.charAt(j)));
			} else {
				alignmentMatrix.get(0).add(0);
			}
		}
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				int[] scoreSet = new int[3];
				scoreSet[0] = alignmentMatrix.get(i).get(j) + scoringMatrix.get(seqX.charAt(i)).get(seqY.charAt(j));
				scoreSet[1] = alignmentMatrix.get(i).get(j + 1) + scoringMatrix.get(seqX.charAt(i)).get('_');
				scoreSet[2] = alignmentMatrix.get(i + 1).get(j) + scoringMatrix.get(seqY.charAt(j)).get('_');
				int maxVal = -(9999);
				for (int k : scoreSet) {
					if (k > maxVal) {
						maxVal = k;
					}
				}
				alignmentMatrix.get(i + 1).add(maxVal);
			}
		}
		return alignmentMatrix;
	}
	
	//Compute global alignment of two sequences
	public void computeGlobalAlignment(String seqX, String seqY, HashMap<Character, HashMap<Character, Integer>> scoringMatrix, 
			ArrayList<ArrayList<Integer>> alignmentMatrix) {
		int iVal = seqX.length();
		int jVal = seqY.length();
		String xPrime = "";
		String yPrime = "";
		int score = alignmentMatrix.get(iVal).get(jVal);
		
		while ((iVal != 0) && (jVal != 0)) {
			int diagonal = alignmentMatrix.get(iVal - 1).get(jVal - 1) + scoringMatrix.get(seqX.charAt(iVal - 1)).get(jVal - 1);
			int vertical = alignmentMatrix.get(iVal - 1).get(jVal) + scoringMatrix.get(seqX.charAt(iVal - 1)).get('_');
			if (alignmentMatrix.get(iVal).get(jVal) == diagonal) {
				xPrime = seqX.charAt(iVal - 1) + xPrime;
				yPrime = seqY.charAt(jVal - 1) + yPrime;
				iVal--;
				jVal--;
			} else {
				if (alignmentMatrix.get(iVal).get(jVal) == vertical) {
					xPrime = seqX.charAt(iVal - 1) + xPrime;
					yPrime = '_' + yPrime;
					iVal--;
				} else {
					xPrime = '_' + xPrime;
					yPrime = seqY.charAt(jVal - 1) + yPrime;
					jVal--;
				}
			}
		}
		while (iVal != 0) {
			xPrime = seqX.charAt(iVal - 1) + xPrime;
			yPrime = '_' + yPrime;
			iVal--;
		}
		while (jVal != 0) {
			xPrime = '_' + xPrime;
			yPrime = seqY.charAt(jVal - 1) + yPrime;
			jVal--;
		}
		
		System.out.println(score);
		System.out.println(xPrime);
		System.out.println(yPrime);
	}
	
	public void computeLocalAlignment(String seqX, String seqY, HashMap<Character, HashMap<Character, Integer>> scoringMatrix, 
			ArrayList<ArrayList<Integer>> alignmentMatrix) {
		int iVal = seqX.length();
		int jVal = seqY.length();
		String xPrime = "";
		String yPrime = "";
		int maxScore = 0;
		
		for (int x = 0; x < seqX.length() + 1; x++) {
			for (int y = 0; y < seqY.length() + 1; y++) {
				int score = alignmentMatrix.get(x).get(y);
				if (maxScore < score) {
					maxScore = score;
					iVal = x;
					jVal = y;
				}
			}
		}
		while (alignmentMatrix.get(iVal).get(jVal) != 0) {
			int diagonal = alignmentMatrix.get(iVal - 1).get(jVal - 1) + scoringMatrix.get(seqX.charAt(iVal - 1)).get(seqY.charAt(jVal - 1));
			int vertical = alignmentMatrix.get(iVal - 1).get(jVal) + scoringMatrix.get(seqX.charAt(iVal - 1)).get('_');
			if (alignmentMatrix.get(iVal).get(jVal) == diagonal) {
				xPrime = seqX.charAt(iVal - 1) + xPrime;
				yPrime = seqY.charAt(jVal - 1) + yPrime;
				iVal--;
				jVal--;
			} else {
				if (alignmentMatrix.get(iVal).get(jVal) == vertical) {
					xPrime = seqX.charAt(iVal - 1) + xPrime;
					yPrime = '_' + yPrime;
					iVal--;
				} else {
					xPrime = '_' + xPrime;
					yPrime = seqY.charAt(jVal - 1) + yPrime;
					jVal--;
				}
			}
		}
		int score = 0;
		for (int i = 0; i < xPrime.length(); i++) {
			score += scoringMatrix.get(xPrime.charAt(i)).get(yPrime.charAt(i));
		}
		
		System.out.println(score);
		System.out.println(xPrime);
		System.out.println(yPrime);
	}
}
