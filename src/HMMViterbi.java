// Emily Behrendt and Katie McCorkell
// eabehr, 1128821 // kmccork, 0822555
// CSE 427 Homework 3
// Thursday, February 20, 2013
// Viterbi Algorithm


import java.io.*;
import java.util.*;

public class HMMViterbi {

	public static final String actg = "ACGT"; // gene index
	public static final String die = "123456";
	
	public static double[][] genTransitions = new double[][] {{.9999, .0001}, {.9999, .0001}, {.01, .99}};
	public static double[] genEState1 = new double[] {.25, .25, .25, .25};
	public static double[] genEState2 = new double[] {.20, .30, .30, .20};

	public static char[] viPath;
	public static char[] genome;
	
	public static final char replacement = 'T'; // if input is not one of these it is converted to a 'T'

	// variables for Viterbi on loaded/fair dice rolls (requires minor changes to code to run on this data)
	public static final String DICE_SEQ = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";
	public static double[][] dieTransitions = new double[][] {{.52, .48}, {.60, .40}, {.17, .83}};
	public static double[][] bigDieTransitions = new double[][] {{.1, .9}, {.9, .1}, {.05, .95}};
	public static double[] dieELoaded = new double[] {.10, .10, .10, .10, .10, .5};
	public static final double fair = 1.0/6.0;
	public static double[] dieEFair = new double[] {fair, fair, fair, fair, fair, fair};
	
	public static void main(String[] args) throws FileNotFoundException {
		readGene();
		cleanGene();

		printStart();
		
		double begin = System.nanoTime();
		double end = 0.0;

		boolean printHits;
		for(int i = 1; i <= 10; i++) {
			printHits = (i == 1 || i == 10);
			System.out.println("Iteration " + i);
			hmmViterbi(genome, genTransitions, genEState1, genEState2);
			processPath(printHits);
			System.out.println("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
			System.out.println();
			if(i == 9) {
				end = System.nanoTime();  
			}
		}
		System.out.println("Time: " + (end - begin) + " nanoseconds");
	}
	
	// start - start of section, end - end of section (inclusive), state2 - if true, state1, if false, state1
	public static void generateNewProbabilities(int start, int end, boolean state2) {
		for(int i = start; i <= end; i++) {
			if(state2) {
				genEState2[actg.indexOf(genome[i])]++;
			} else {
				genEState1[actg.indexOf(genome[i])]++;
			}
		}
	}
	
	// This method processes a given Viterbi path.
	// It identifies the hits (contiguous sequences of State 2) and if the boolean parameter is true, prints them
	// It also does the E-step of the EM algorithm, aka:
	// For the entire path, it recalculates the emission and transition probabilities,
	// and replaces them in the global variables 
	// Note: It doesn't change the begin state probabilities
	public static void processPath(boolean printHits) {
		// reset probabilities
		genEState1 = new double[]{0.0, 0.0, 0.0, 0.0};
		genEState2 = new double[]{0.0, 0.0, 0.0, 0.0};
		genTransitions = new double[][]{{.9999, .0001}, {0.0, 0.0}, {0.0, 0.0}};

		// state2 = 1, state1 = 0, looking for continuous sequences of 1s	  
		System.out.println("Lengths and Locations of All Hits (start and end both inclusive)");
		// whether currently in a sequence of 1s
		boolean inSeq = (viPath[0] == '1');	  
		int numHits = 0;
		int start, end, length;
		for(int i = 0; i < viPath.length; i++) {
			start = i;
			while(i < viPath.length && viPath[i] == '0') {
				if(inSeq) {
					// increment state2->state1
					genTransitions[2][0]++;
				} else {
					// increment state1->state1
					genTransitions[1][0]++;
				}
				inSeq = false;
				i++;
			}

			end = i-1;
			generateNewProbabilities(start, end, false);
			start = i;
			while(i < viPath.length && viPath[i] == '1') {
				if(!inSeq) {
					// increment state1->state2
					genTransitions[1][1]++;
				} else {
					// increment state2->state2
					genTransitions[2][1]++;
				}
				inSeq = true;
				i++;
			}
			end = i-1;
			length = end - start + 1;
			if(inSeq) {
				numHits++;
				if(printHits) {
					System.out.println("Start - End: [" + (start+1) + ", " + (end+1) + "]\t" + "Length: " + length);
				}
				generateNewProbabilities(start, end, true);
			}
		}
		System.out.println("Number of hits: " + numHits + "\n");

		double state1totaltrans = genTransitions[1][0] + genTransitions[1][1];  
		double state2totaltrans = genTransitions[2][0] + genTransitions[2][1];

		genTransitions[1][0] = genTransitions[1][0] / state1totaltrans;
		genTransitions[1][1] = genTransitions[1][1] / state1totaltrans;
		genTransitions[2][0] = genTransitions[2][0] / state2totaltrans;
		genTransitions[2][1] = genTransitions[2][1] / state2totaltrans;

		System.out.println("Transition Probabilities");
		System.out.println("- - - - - - - - - - - - -");
		System.out.println("Transitions\t\tState1\t\tState2");
		System.out.println("State1\t\t" + genTransitions[1][0] + "\t" + genTransitions[1][1]);
		System.out.println("State2\t\t" + genTransitions[2][0] + "\t" + genTransitions[2][1]);
		System.out.println();

		double state1total = 0;
		double state2total = 0;
		for(int i = 0; i < 4; i++) {
			state1total += genEState1[i];
			state2total += genEState2[i];
		}
		for(int i = 0; i < 4; i++) {
			genEState1[i] = genEState1[i]/state1total;
			genEState2[i] = genEState2[i]/state2total;
		}

		System.out.println("\t\t\tA\t\t\tC\t\t\tG\t\t\tT");
		System.out.println("State1 Probabilities: " + Arrays.toString(genEState1));
		System.out.println("State2 Probabilities: " + Arrays.toString(genEState2));
	}

	// trans = transition = "a" 
	// e = emitL = emissions
	// Calculates the overall Viterbi path for the sequence
	// Stores this path in global variable viPath
	public static void hmmViterbi(char[] input, double[][] trans, double[] emitL, double[] emitF){
		// this output grid looks like 
		// output[0][i] 11  --> probability of emission from state1 given previously state1
		// output[1][i] 21  --> probability of emission from state1 given previously state2
		// output[2][i] Max of State1 state possibilities 
		// output[3][i] 12  --> probability of emission from state2 given previously state1
		// output[4][i] 22  --> probability of emission from state2 given previously state2
		// output[5][i] Max of State2 state possibilities  
		double[][] output = new double[6][input.length];
		for(int r = 0; r < 3; r++) {
			// begin -> state1 transition
			output[r][0] = Math.log(trans[0][0]) + Math.log(emitL[actg.indexOf(input[0])]);
		} 
		for(int r = 3; r < 6; r++) {
			// begin -> state2 transition
			output[r][0] = Math.log(trans[0][1]) + Math.log(emitF[actg.indexOf(input[0])]);
		} 

		int[][] path = new int[2][input.length];
		for(int i = 1; i < input.length; i++) {
			// s1->s1 transition
			output[0][i] = (output[2][i-1]) + Math.log(trans[1][0]) + Math.log(emitL[actg.indexOf(input[i])]);
			// s2->s1 transition
			output[1][i] = (output[5][i-1]) + Math.log(trans[2][0]) + Math.log(emitL[actg.indexOf(input[i])]);
			output[2][i] = output[0][i];
			if(output[1][i] > output[0][i]) {
				path[0][i] = 1;
				output[2][i] = output[1][i];
			}
			// s1->s2 transition
			output[3][i] = (output[2][i-1]) + Math.log(trans[1][1]) + Math.log(emitF[actg.indexOf(input[i])]);
			// s2->s2 transition
			output[4][i] = (output[5][i-1]) + Math.log(trans[2][1]) + Math.log(emitF[actg.indexOf(input[i])]);
			output[5][i] = output[3][i];
			if(output[4][i] > output[3][i]) {
				path[1][i] = 1;
				output[5][i] = output[4][i];
			}
		}

		double viterbiPathProb = Math.max(output[5][input.length-1], output[2][input.length-1]);

		System.out.println("Log probability of the overall Viterbi path: " + viterbiPathProb + "\n");

		viPath = new char[input.length];

		// traceback to Get the viterbi path 
		int i = input.length -1;
		int c;
		if(output[5][i] > output[2][i]){
			viPath[i] = '1';
			c = path[1][i]; 
		} else {
			viPath[i] = '0';
			c = path[0][i];
		} 
		for(i = input.length - 2; i >= 0; i--) {
			if(c == 0) {
				viPath[i] = '0';
			} else {
				viPath[i] = '1';
			}
			c = path[c][i];
		}
	} 

	// Replaces non valid nucleotides with the replacement nucleotide 
	private static void cleanGene() {
		for (int i = 0; i < genome.length; i++) {
			if (genome[i] != 'A' && genome[i] != 'C' && genome[i] != 'G'
					&& genome[i] != 'T') {
				genome[i] = 'T';
			}
		}
	}

	// Reads the gene file into a char array and stores it in global variable  
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
			sequence += line;
		}
		genome = sequence.toCharArray();
	}
	
	public static void printStart() {
		System.out.println("Starting Transition Probabilities");
		System.out.println("- - - - - - - - - - - - - - - - -");
		System.out.println("Transitions\tState1\t\tState2");
		System.out.println("Begin\t\t" + genTransitions[0][0] + "\t\t" + genTransitions[0][1]);
		System.out.println("State1\t\t" + genTransitions[1][0] + "\t\t" + genTransitions[1][1]);
		System.out.println("State2\t\t" + genTransitions[2][0] + "\t\t" + genTransitions[2][1]);
		System.out.println();
	}
}
