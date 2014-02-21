import java.io.*;
import java.util.*;

public class HMMViterbi {
	
  public static final char replacement = 'T'; // if input is not one of these it is converted to a 'T'
  public static Set<String> stopCodons;

  public static final String DICE_SEQ = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";

  public static final String actg = "ACGT"; // gene index
  public static final String die = "123456";
  
  public static double[][] genTransitions = new double[][] {{.9999, .0001}, {.9999, .0001}, {.01, .99}};
  public static double[][] dieTransitions = new double[][] {{.52, .48}, {.60, .40}, {.17, .83}};
  public static double[][] bigDieTransitions = new double[][] {{.1, .9}, {.9, .1}, {.05, .95}};
  
  public static double[] dieELoaded = new double[] {.10, .10, .10, .10, .10, .5};
  public static final double fair = 1.0/6.0;
  public static double[] dieEFair = new double[] {fair, fair, fair, fair, fair, fair};
  
  public static double[] genEState1 = new double[] {.25, .25, .25, .25};
  public static double[] genEState2 = new double[] {.20, .30, .30, .20};

  public static char[] viPath;
  public static char[] genome;
  
  public static void main(String[] args) throws FileNotFoundException {

    genome = readGene();
    cleanGene(genome);
// remove genome as parameter!
    hmmViterbi(genome, genTransitions, genEState1, genEState2);
    
    processPath();
    
    
    
   // char[] diceSeq = prepDiceSeq();
   // char[] shortDie = "666666".toCharArray();
   // hmmViterbi(diceSeq, bigDieTransitions, dieELoaded, dieEFair);
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
  
  public static void processPath() {
	  // reset probabilities
	  genEState1 = new double[]{0.0, 0.0, 0.0, 0.0};
	  genEState2 = new double[]{0.0, 0.0, 0.0, 0.0};
	  
	  //public static double[][] genTransitions = new double[][] {{.9999, .0001}, {.9999, .0001}, {.01, .99}};
	  genTransitions = new double[][]{{.9999, .0001}, {0.0, 0.0}, {0.0, 0.0}};
	  //double 
	  
	  // state2 = 1, state1 = 0, looking for continuous sequences of 1s	  
	  System.out.println("Lengths and Locations of All Hits");
	  boolean inSeq = (viPath[0] == '1'); // whether currently in a sequence of 1s	  
	  int numHits = 0;
	  int start, end, length;
	  for(int i = 0; i < viPath.length; i++) {
		  start = i;
		  while(i < viPath.length && viPath[i] == '0') {
			  if(inSeq) {
				  //increment state 2 ->state1
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
			  System.out.println("Start: " + (start+1) + "\t" + "End: " + (end+1) + "\t" + "Length: " + length);
			  generateNewProbabilities(start, end, true);
		  }
	  }
	  System.out.println("Number of hits: " + numHits);
	  
	  		// maybe same as state1 total, state 2 total???
	  double state1totaltrans = genTransitions[1][0] + genTransitions[1][1];  
	  double state2totaltrans = genTransitions[2][0] + genTransitions[2][1];
	  
	  genTransitions[1][0] = genTransitions[1][0] / state1totaltrans;
	  genTransitions[1][1] = genTransitions[1][1] / state1totaltrans;
	  genTransitions[2][0] = genTransitions[2][0] / state2totaltrans;
	  genTransitions[2][1] = genTransitions[2][1] / state2totaltrans;
	  
	  print2Array(genTransitions, 3, 2);
	  
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
	  
	  System.out.println("State1 Probabilities: " + Arrays.toString(genEState1));
	  System.out.println("State2 Probabilities: " + Arrays.toString(genEState2));
  }
  
  //trans = transition = "a" 
  //e = emitL = emissions
  public static void hmmViterbi(char[] input, double[][] trans, double[] emitL, double[] emitF){
    // this output grid looks like 
    // output[0][i] LL  --> probability of role in loaded state given previous state was loaded die
    // output[1][i] FL  --> probability of role in loaded state given previous state was fair die
    // output[2][i] Max of Loaded state possibilities 
    // output[3][i] LF  --> probability of role in fair state given previous state was loaded die 
    // output[4][i] FF  --> probability of role in fair state given previous state was fair die
    // output[5][i] Max of Fair state possibilities  
    double[][] output = new double[6][input.length];
    for(int r=0; r<3; r++) {
      output[r][0] = Math.log(trans[0][0]) + Math.log(emitL[actg.indexOf(input[0])]); //B --> L transition
    } 
    for(int r=3; r<6; r++) {
      output[r][0] = Math.log(trans[0][1]) + Math.log(emitF[actg.indexOf(input[0])]); //B --> F transition
    } 

    int[][] path = new int[2][input.length];
    for(int i=1; i<input.length; i++) {
      output[0][i] = (output[2][i-1]) + Math.log(trans[1][0]) + Math.log(emitL[actg.indexOf(input[i])]); //LL transition --> Loaded die to loaded die
      output[1][i] = (output[5][i-1]) + Math.log(trans[2][0]) + Math.log(emitL[actg.indexOf(input[i])]); // FL transition
      output[2][i] = output[0][i];
      if(output[1][i] > output[0][i]) {
        path[0][i] = 1;
        output[2][i] = output[1][i];
      }

      output[3][i] = (output[2][i-1]) + Math.log(trans[1][1]) + Math.log(emitF[actg.indexOf(input[i])]); //LF transition
      output[4][i] = (output[5][i-1]) + Math.log(trans[2][1]) + Math.log(emitF[actg.indexOf(input[i])]); // FF transition
      output[5][i] = output[3][i];
      if(output[4][i] > output[3][i]) {
        path[1][i] = 1;
        output[5][i] = output[4][i];
      }
    }
    
    double viterbiPathProb = Math.max(output[5][input.length-1], output[2][input.length-1]);
    
    System.out.println("Log probability of the overall Viterbi path: " + viterbiPathProb + "\n");
  
    //print2Array(output, 6, input.length);

    viPath = new char[input.length];

    // traceback
    // To Get the viterbi path 
    int i = input.length -1;
    int c;
    if(output[5][i] > output[2][i]){
      viPath[i] = '1';
      c = path[1][i]; 
    } else {
      viPath[i] = '0';
      c = path[0][i];
    } 
    for( i=input.length-2; i>=0; i--) {
     if(c == 0) {
        viPath[i] = '0';
      } else {
        viPath[i] = '1';
      }
      c = path[c][i];
    }
  } 

  // Preps the sequence for testing the casino model 
  public static char[] prepDiceSeq() {
    char[] diceSeq = DICE_SEQ.toCharArray();
    return diceSeq;
  }

  public static double round(double value) {
    return (double)Math.round(value * 10000) / 10000;
  }
  // Prints a 2-d array pretty decently 
  public static void print2Array(double[][] array, int rows, int columns){
    for(int i=0; i<rows; i++) {
      for(int j=0; j<columns; j++) {
        System.out.print((array[i][j]) + "\t");
      }
      System.out.println();
    }
  }

	// Replaces non valid nucleotides with the replacement nucleotide 
	private static char[] cleanGene(char[] genome) {
		for(int i = 0; i<genome.length; i++) {
			String gi = "" + genome[i]; //contains method takes a string, not a char
			if(!actg.contains(gi)) {
				genome[i] = replacement;
			}
		}
		return genome;
	}

  // Reads the gene file into a char array and returns it 
	public static char[] readGene() throws FileNotFoundException{
	    File f = new File("FNA");
	    Scanner s = new Scanner(f);
	    char[] genome = new char[0];
	    String longLine = "";
	    while (s.hasNextLine()) {
	      Scanner lineScanner = new Scanner(s.nextLine());
	      while(lineScanner.hasNext()) {
	    	  String line = lineScanner.nextLine();
          line = line.toUpperCase(); 
          if (line.charAt(0) != '>') {
            line = line.trim();
            longLine = longLine + line; 
          }
	      }
	    }
		genome = longLine.toCharArray();
	  return genome;
	}
}
