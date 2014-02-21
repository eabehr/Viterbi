import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class BaumWelch {

	public static double[][] transitionProbs = new double[][] {{.9999, .0001}, {.9999, .0001}, {.01, .99}};
	public double[][] emissionProbs = new double[][] {{.25, .25, .25, .25}, {.20, .30, .30, .20}};
	
	public static char[] genome;
	
	public static void main(String[] args) throws FileNotFoundException {
		readGene();
		cleanGene();
		
		double[][] forward = forward();
		double[][] backward;
		
		for(int i = 0; i < 10; i++) {
			
		}
	}

	public static double[][] forward() {
		double[][] forward = new double[2][genome.length];
		
		return forward;
	}
	
	// Replaces non valid nucleotides with the replacement nucleotide 
	private static void cleanGene() {
		for(int i = 0; i < genome.length; i++) {
			if(genome[i] != 'A' && genome[i] != 'C' && genome[i] != 'G' && genome[i] != 'T') {
				genome[i] = 'T';
			}
		}
	}

	// Reads the gene file into a char array and returns it 
	public static void readGene() throws FileNotFoundException{
		File f = new File("FNA");
		Scanner s = new Scanner(f);
		genome = new char[0];
		String sequence = "";
		// first line of FNA file not part of sequence
		String line = s.nextLine();
		while (s.hasNextLine()) {
			line = s.nextLine();	
			line = line.toUpperCase(); 
			line = line.trim();
			sequence += line;; 
		}
		genome = sequence.toCharArray();
	}
}
