import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

public class BaumWelch {

	public static String actg = "ACGT"; // gene index
	public static final String die = "123456";

	public static final String DICE_SEQ = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";

	public static final char[] dieRolls = DICE_SEQ.toCharArray(); 
	public static double[][] transitionProbs = new double[][] {
			{ .9999, .0001 }, { .9999, .0001 }, { .01, .99 } };
	public static double[][] emissionProbs = new double[][] {
			{ .25, .25, .25, .25 }, { .20, .30, .30, .20 } };
	
	public static double[][] dieTransitionProbs = new double[][] {{.1, .9}, {.9, .1}, {.05, .95}};
	public static final double fair = 1.0/6.0;
	public static double[][] dieEmissionProbs = new double[][] {{.10, .10, .10, .10, .10, .5}, {fair, fair, fair, fair, fair, fair}};

	public static double pathProb;

	public static char[] input; 
	public static char[] genome;
	public static char[] viPath; 
	public static double[][] backward; 
	public static double[][] forward; 
	public static double[][] forwardCompute; 
	public static double[][] originalTransitions; 
	public static double[][] originalEmissions; 


	public static void main(String[] args) throws FileNotFoundException {
		readGene();
		cleanGene();
		
		printStart();
		double begin = System.nanoTime();
		double end = 0.0;
		input = genome; 		
		boolean diceExample = false; 
		if (diceExample){
			input = DICE_SEQ.toCharArray();
			genome = input; 
			transitionProbs = dieTransitionProbs;
			emissionProbs = dieEmissionProbs; 
			actg = die; 
		}
		originalEmissions = emissionProbs; 
		originalTransitions = transitionProbs;
		boolean printHits;
		for (int i = 1; i <= 10; i++) {
			printHits = (i == 1 || i == 10);
			System.out.println("Iteration " + i);
			if(diceExample) {
				hmmViterbi(dieRolls, dieTransitionProbs, dieEmissionProbs[0], dieEmissionProbs[1], die); 
			}else {
				hmmViterbi(genome, transitionProbs, emissionProbs[0], emissionProbs[1], actg); 
			}
			forwardAlgorithm();
			backwardAlgorithm();
			processPath(printHits, diceExample);
			System.out
					.println("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
			System.out.println();
			if (i == 9) {
				end = System.nanoTime();
			}
			originalTransitions = transitionProbs;
			originalEmissions = emissionProbs; 
		}
		System.out.println("Time: " + (end - begin) + " nanoseconds");

	}

	// The number of transitions from k to l
	// is the sum of the probability of transition from k to l
	// for each nucleotide (that is controlled by different method from this
	// one)
	private static void updateTransitionProbs(int i, int seqLength) {
		for (int row = 1; row < 3; row++) {
			for (int column = 0; column < 2; column++) {
				transitionProbs[row][column] += probabilityOfTransitionKL(
						row, column, i, seqLength);
			}
		}

	}

	// k is the state you're in either 1 or 2
	// L is the state you're transitioning to either 1 or 2
	private static double probabilityOfTransitionKL(int k, int l, int i, int seqLength) {
		if (i >= seqLength -1){
			return 0;
		}
		int xiplusOneIndex = actg.indexOf(input[i + 1]);
		double f = forward[k-1][i];
		double t = originalTransitions[k][l];
		double e = originalEmissions[l][xiplusOneIndex];
		double b = backward[l][i+1];
		double numerator = f * t * e * b;
		double p = numerator / pXGivenTheta();
		return p;
	}

	private static double pXGivenTheta() {
		return pathProb;
	}

	public static void printStart() {
		System.out.println("Starting Transition Probabilities");
		System.out.println("- - - - - - - - - - - - - - - - -");
		System.out.println("Transitions\tState1\t\tState2");
		System.out.println("Begin\t\t" + transitionProbs[0][0] + "\t\t"
				+ transitionProbs[0][1]);
		System.out.println("State1\t\t" + transitionProbs[1][0] + "\t\t"
				+ transitionProbs[1][1]);
		System.out.println("State2\t\t" + transitionProbs[2][0] + "\t\t"
				+ transitionProbs[2][1]);
		System.out.println();
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

	// Reads the gene file into a char array and returns it
	public static void readGene() throws FileNotFoundException {
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

	/*
	 * Given two probabilities x and y, represented by their logs lx, ly, return
	 * the log of their sum log(x+y) = log(exp(lx) + exp(ly)).
	 * 
	 * Assume log(0) is represented by NaN.
	 * 
	 * The "lx > ly" trick is some protection from underflow: log(a+b) =
	 * log(a(1+b/a)) = log(a)+log(1+b/a), which will be most accurate when b/a <
	 * 1.
	 */
	public static double log_of_sum_of_logs(double lx, double ly) {
		if (isnan(lx))
			return ly;
		if (isnan(ly))
			return lx;
		if (lx > ly) {
			return lx + Math.log(1 + Math.exp(ly - lx));
		} else {
			return ly + Math.log(1 + Math.exp(lx - ly));
		}
	}

	public static boolean isnan(Double x) {
		return x.isNaN();
	}

	// This method processes a given posterior path.
	// It identifies the hits and if the boolean parameter is true, prints them
	// A "hit" is a (contiguous) subsequence with state 2 posterior probability
	// above one half
	// -----
	// It updates the transition probabilities using the Baum-Welch approach
	// It recalculates the emission probabilities in the same way as the Viterbi
	// approach
	// and replaces both in the global variables
	// Note: It doesn't change the begin state probabilities
	public static void processPath(boolean printHits, boolean diceExample) {
		
		int seqLength = input.length; 

		// reset probabilities
		emissionProbs = new double[][] { { 0.0, 0.0, 0.0, 0.0},
				{ 0.0, 0.0, 0.0, 0.0} };
		transitionProbs = new double[][] { { .9999, .0001 }, { 0.0, 0.0 },
				{ 0.0, 0.0 } };

		if(diceExample) {
			emissionProbs = new double[][] { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
		}
		// state2 = 1, state1 = 0, looking for continuous sequences of 1s
		System.out
				.println("Lengths and Locations of All Hits (start and end both inclusive)");
		// whether currently in a sequence of 1s

		boolean inSeq = (viPath[0] == '1');
		int numHits = 0;
		int start, end, length;
		for (int i = 0; i < input.length; i++) {
			start = i;
			while (i < input.length && viPath[i] == '0') {
				updateTransitionProbs(i, seqLength);
				inSeq = false;
				i++;
			}

			end = i - 1;
			generateNewProbabilities(start, end, false);
			start = i;
			while (i < input.length && viPath[i] == '1') {
				updateTransitionProbs(i, seqLength);
				inSeq = true;
				i++;
			}
			end = i - 1;
			length = end - start + 1;
			if (inSeq) {
				numHits++;
				if (printHits) {
					System.out.println("Start - End: [" + (start + 1) + ", "
							+ (end + 1) + "]\t" + "Length: " + length);
				}
				generateNewProbabilities(start, end, true);
			}
		}
		System.out.println("Number of hits: " + numHits + "\n");

		double state1totaltrans = transitionProbs[1][0] + transitionProbs[1][1];
		double state2totaltrans = transitionProbs[2][0] + transitionProbs[2][1];

		transitionProbs[1][0] = transitionProbs[1][0] / state1totaltrans;
		transitionProbs[1][1] = transitionProbs[1][1] / state1totaltrans;
		transitionProbs[2][0] = transitionProbs[2][0] / state2totaltrans;
		transitionProbs[2][1] = transitionProbs[2][1] / state2totaltrans;

		System.out.println("Transition Probabilities");
		System.out.println("- - - - - - - - - - - - -");
		System.out.println("Transitions\t\tState1\t\tState2");
		System.out.println("State1\t\t" + transitionProbs[1][0] + "\t"
				+ transitionProbs[1][1]);
		System.out.println("State2\t\t" + transitionProbs[2][0] + "\t"
				+ transitionProbs[2][1]);
		System.out.println();

		double state1total = 0;
		double state2total = 0;
		for (int i = 0; i < emissionProbs[0].length; i++) {
			state1total += emissionProbs[0][i];
			state2total += emissionProbs[1][i];
		}
		for (int i = 0; i < emissionProbs[0].length; i++) {
			emissionProbs[0][i] = emissionProbs[0][i] / state1total;
			emissionProbs[1][i] = emissionProbs[1][i] / state2total;
		}

		System.out.println("\t\t\tA\t\t\tC\t\t\tG\t\t\tT");
		System.out.println("State1 Probabilities: "
				+ Arrays.toString(emissionProbs[0]));
		System.out.println("State2 Probabilities: "
				+ Arrays.toString(emissionProbs[1]));
	}

	// start - start of section, end - end of section (inclusive), state2 - if
	// true, state1, if false, state1
	public static void generateNewProbabilities(int start, int end,
			boolean state2) {
		for (int i = start; i <= end; i++) {
			if (state2) {
				emissionProbs[1][actg.indexOf(genome[i])]++;
			} else {
				emissionProbs[0][actg.indexOf(genome[i])]++;
			}
		}
	}
	
	public static void forwardAlgorithm(){
		// LL  --> probability of role in loaded state given previous state was loaded die
		// FL  --> probability of role in loaded state given previous state was fair die
		// LF  --> probability of role in fair state given previous state was loaded die 
		// FF  --> probability of role in fair state given previous state was fair die
		forward = new double[2][input.length];
		forwardCompute = new double[4][input.length]; 
		// Begin State Transition Probabilities 
		for(int r = 0; r < 2; r++) {
			forward[r][0] = Math.log(transitionProbs[0][r]) + Math.log(emissionProbs[r][actg.indexOf(input[0])]); 
		} 

		for(int i = 1; i < input.length; i++) {
			forwardCompute[0][i] = (forwardCompute[0][i-1]) + Math.log(transitionProbs[1][0]) + Math.log(emissionProbs[0][actg.indexOf(input[i])]); //LL transition --> Loaded die to loaded die
			forwardCompute[1][i] = (forwardCompute[1][i]) + Math.log(transitionProbs[2][0]) + Math.log(emissionProbs[0][actg.indexOf(input[i])]); // FL transition
			forward[0][i] = log_of_sum_of_logs(forwardCompute[0][i], forwardCompute[1][i]); 

			forwardCompute[2][i] = (forwardCompute[2][i-1]) + Math.log(transitionProbs[1][1]) + Math.log(emissionProbs[1][actg.indexOf(input[i])]); //LF transition
			forwardCompute[3][i] = (forwardCompute[3][i]) + Math.log(transitionProbs[2][1]) + Math.log(emissionProbs[1][actg.indexOf(input[i])]); // FF transition
			forward[1][i] = log_of_sum_of_logs(forwardCompute[2][i], forwardCompute[3][i]);
		}
	} 
	
	
	public static void backwardAlgorithm(){
		
		// LL  --> probability of role in loaded state given previous state was loaded die
		// FL  --> probability of role in loaded state given previous state was fair die
		// LF  --> probability of role in fair state given previous state was loaded die 
		// FF  --> probability of role in fair state given previous state was fair die
		backward = new double[2][input.length];


		for(int i = input.length-2; i > 0; i--) {
			int xIndex = actg.indexOf(input[i+1]); 

			double LL = (backward[0][i+1]) + Math.log(transitionProbs[1][0]) + Math.log(emissionProbs[0][xIndex]); 
			double LF = (backward[1][i+1]) + Math.log(transitionProbs[1][1]) + Math.log(emissionProbs[1][xIndex]);
			//0 is loaded state 
			backward[0][i] = log_of_sum_of_logs(LL, LF); 

			double FL = (backward[0][i+1]) + Math.log(transitionProbs[2][0]) + Math.log(emissionProbs[0][xIndex]);
			double FF = (backward[1][i+1]) + Math.log(transitionProbs[2][1]) + Math.log(emissionProbs[1][xIndex]); 
			backward[1][i] = log_of_sum_of_logs(FL, FF);
		}
	} 
	
	//trans = transition = "a" 
		//e = emitL = emissions
		// Calculates the overall Viterbi path for the sequence
		// Stores this path in global variable viPath
		public static void hmmViterbi(char[] input, double[][] trans, double[] emitL, double[] emitF, String code){
			// this output grid looks like 
			// output[0][i] LL  --> probability of role in loaded state given previous state was loaded die
			// output[1][i] FL  --> probability of role in loaded state given previous state was fair die
			// output[2][i] Max of Loaded state possibilities 
			// output[3][i] LF  --> probability of role in fair state given previous state was loaded die 
			// output[4][i] FF  --> probability of role in fair state given previous state was fair die
			// output[5][i] Max of Fair state possibilities  
			double[][] output = new double[6][input.length];
			for(int r = 0; r < 3; r++) {
				output[r][0] = Math.log(trans[0][0]) + Math.log(emitL[code.indexOf(input[0])]); //B --> L transition
			} 
			for(int r = 3; r < 6; r++) {
				output[r][0] = Math.log(trans[0][1]) + Math.log(emitF[code.indexOf(input[0])]); //B --> F transition
			} 

			int[][] path = new int[2][input.length];
			for(int i = 1; i < input.length; i++) {
				int xIndex = code.indexOf(input[i]);
				output[0][i] = (output[2][i-1]) + Math.log(trans[1][0]) + Math.log(emitL[xIndex]); //LL transition --> Loaded die to loaded die
				output[1][i] = (output[5][i-1]) + Math.log(trans[2][0]) + Math.log(emitL[xIndex]); // FL transition
				output[2][i] = output[0][i];
				if(output[1][i] > output[0][i]) {
					path[0][i] = 1;
					output[2][i] = output[1][i];
				}

				output[3][i] = (output[2][i-1]) + Math.log(trans[1][1]) + Math.log(emitF[code.indexOf(input[i])]); //LF transition
				output[4][i] = (output[5][i-1]) + Math.log(trans[2][1]) + Math.log(emitF[code.indexOf(input[i])]); // FF transition
				output[5][i] = output[3][i];
				if(output[4][i] > output[3][i]) {
					path[1][i] = 1;
					output[5][i] = output[4][i];
				}
			}

			pathProb = Math.max(output[5][input.length-1], output[2][input.length-1]);

			System.out.println("Log probability of the overall Viterbi path: " + pathProb + "\n");

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
}
