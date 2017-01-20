import java.util.*;
import java.io.*;

/*
 * A Java implementation of a dynamic programming sequence alignment solution.
 * 
 * This class will read in 2 DNA sequences, a match score, mismatch score, and space score 
 * as command line arguments, and apply the Needleman-Wunsch algorithm to produce the 
 * optimal alignment of the 2 DNA sequences. It will also output the alignment score. A sample input is provided.
 * See: sample_input.txt
 * 
 * Sample output for running this file with command line arguments: sample_input.txt 1 -1 -1
 * 
 * 			CCGCACGTACCAGGCTGTCTCTCGGTTAGAC-AT-AA-ATCTCTTCGCG--AGGAGGCAC-C-GCCCGGTGGTGA-GCGCTTAAGGCAG-GTAGTTGGCCG-AATGACAA                     
 *                    CAGA-AGTCGATGAAGA-C-CTAGGCGATAACAGACACACAGCCCAG-GCT-ATGCGCCGAGGG-AGAATCCTT--CCGCTA-GACGTTGTTGACGGTCTTCAAGTACT
 *			
 *			12
 * 
 * @author John Meredith
 * 
 */
public class AlignSequences {

	private static int m; //match score
	private static int h; //mismatch score
	private static int s; //space score
	private static String i = ""; //row sequence
	private static String j = ""; //column sequence
	private static int maxScoreIndex;
	private static int alignmentScore;

	public static int[][] buildMatrix(){
		int[][] DPmatrix = new int[i.length()+1][j.length()+1]; //create dynamic programming matrix
		
		for(int rows = 1; rows <= i.length(); rows++) //initialize first item in each row to 0
			DPmatrix[rows][0] = 0; 	
		
		for(int cols = 1; cols <= j.length(); cols++) //initialize first item in each column to (col*space score)
			DPmatrix[0][cols] = (cols * s);
		
		//fill in matrix using Needleman-Wunsch
		for(int x = 1; x <= i.length(); x++){
			for(int y = 1; y <= j.length(); y++){
				
				int diag = DPmatrix[x-1][y-1] + (i.charAt(x-1) == j.charAt(y-1) ? m : h);
				int left = DPmatrix[x][y-1] + s;
				int up = DPmatrix[x-1][y] + s;

				DPmatrix[x][y] = Math.max(Math.max(diag, up), left);
			}
		}
		
		maxScoreIndex = findMaxScoreIndex(DPmatrix[i.length()]);
		alignmentScore = DPmatrix[i.length()][maxScoreIndex];
		return DPmatrix; 
	}
	
	/* Trace back through the dynamic programming matrix to obtain the alignment with the maximum score.
	 * This traceback uses the "highroad" alignment for tiebreakers.
	 */
	public static void traceback(int[][] matrix){
		int rowCount = i.length();
		int colCount = maxScoreIndex;
		String seqA = "";
		String seqB = "";
		boolean extraJ = false;
		int extraJCount = 0;
		
		if(colCount != j.length()){
			extraJ = true;
			extraJCount = j.length()-colCount;
		}
		
		while(rowCount > 0 && colCount > 0){
			if(matrix[rowCount][colCount] - s == matrix[rowCount-1][colCount]){
				seqA += i.charAt(rowCount-1);
				seqB += "-";
				rowCount--;

			}else if(matrix[rowCount][colCount] - m == matrix[rowCount-1][colCount-1] && 
					(i.charAt(rowCount-1) == j.charAt(colCount-1))){
				seqA += i.charAt(rowCount-1);
				seqB += j.charAt(colCount-1);
				rowCount--;
				colCount--;

			}else if(matrix[rowCount][colCount] - h == matrix[rowCount-1][colCount-1]){
				seqA += i.charAt(rowCount-1);
				seqB += j.charAt(colCount-1);
				rowCount--;
				colCount--;

			}else if(matrix[rowCount][colCount] - s == matrix[rowCount][colCount-1]){
				seqA += "-";
				seqB += j.charAt(colCount-1);
				colCount--;

			}else{
				System.out.println("Something went wrong, you shouldn't be here.");
			}
		}
		
		if(colCount > 0){
			for(int x = colCount; x > 0; x--){
				seqA += "-";
				seqB += j.charAt(x-1);
			}
		}
		if(rowCount > 0){
			for(int x = rowCount; x > 0; x--){
				seqA += i.charAt(x-1);
				seqB += " ";
			}
		}
		
		seqA = new StringBuilder(seqA).reverse().toString();
		seqB = new StringBuilder(seqB).reverse().toString();
		
		if(extraJ){
			for(int x = extraJCount; x > 0; x--){
				seqB += j.charAt(j.length()-x);
				seqA += " ";
			}
		}
		
		System.out.println(seqA);
		System.out.println(seqB);

	}

	public static int findMaxScoreIndex(int[] a){
		int maxIndex = 0;
		for(int i = 0; i < a.length; i++){
			if(a[maxIndex] <= a[i]){
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	public static void main(String[] args){
		
		//file read
		if(args.length > 0){
			File file = new File(args[0]);
			try{
				Scanner in = new Scanner(file);
				//input file should have 2 lines
				int c = 0;
				while(in.hasNextLine()){
					String line = in.nextLine();
					if(c == 0)
						i += line;
					if(c == 1)
						j += line;
					c++;
				}

			}catch(FileNotFoundException e){
				System.out.println("File not found.");
			}
			m = Integer.parseInt(args[1]);
			h = Integer.parseInt(args[2]);
			s = Integer.parseInt(args[3]);
		}
		
		traceback(buildMatrix());
		System.out.println(alignmentScore);
	}
}
