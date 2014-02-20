import java.io.*;
import java.util.*;

public class HMMViterbi {
	
	public static final char replacement = 'T'; // if input is not one of these it is converted to a 'T'
	public static Set<String> stopCodons;

  public static final String DICE_SEQ = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";

	public static final String actg = "ACTG"; // gene index
  public static final String die = "123456";
  
  public static double[][] genTransitions = new double[][] {{.9999, .0001},
                                                            {.9999, .0001},
                                                            {.01, .99}};

  public static double[][] dieTransitions = new double[][] {{.52, .48},
                                                            {.60, .40}, 
                                                            {.17, .83}};

  public static double[][] bigDieTransitions = new double[][] {{.1, .9},
                                                            {.9, .1}, 
                                                            {.05, .95}};

  public static double[] genEState1 = new double[] {.25, .25, .25, .25};
  public static double[] genEState2 = new double[] {.20, .30, .30, .20};

  public static double[] dieELoaded = new double[] {.10, .10, .10, .10, .10, .5};
  public static final double fair = 1.0/6.0;
  public static double[] dieEFair = new double[] {fair, fair, fair, fair, fair, fair};

	public static void main(String[] args) throws FileNotFoundException {
/*
    char[] genome = readGene();
    genome = cleanGene(genome);
*/
    
    char[] diceSeq = prepDiceSeq();
    char[] shortDie = "666666".toCharArray();
    hmmViterbi(diceSeq, bigDieTransitions, dieELoaded, dieEFair);
  }
  
  //trans = transition = "a" 
  //e = emitL = emissions
  public static void hmmViterbi(char[] input, double[][] trans, double[] emitL, double[] emitF){
    // this output grid looks like 
    // output[0][i] LL  --> probability of role in loaded state given previous state was loaded die
    // output[1][i] FL  --> probability of role in loaded state given previous state was fair die
    // output[2][i] Max of Loaded state posibilities 
    // output[3][i] LF  --> probability of role in fair state given previous state was loaded die 
    // output[4][i] FF  --> probability of role in fair state given previous state was fair die
    // output[5][i] Max of Fair state possibilities  
    double[][] output = new double[6][input.length];
    for(int r=0; r<3; r++) {
      output[r][0] = (trans[0][0]) * (emitL[die.indexOf(input[0])]); //B --> L transition
    } 
    for(int r=3; r<6; r++) {
      output[r][0] = (trans[0][1]) * (emitF[die.indexOf(input[0])]);// B--> F transition
    } 

    int[][] path = new int[2][input.length];
    for(int i=1; i<output.length; i++) {
      output[0][i] = (output[2][i-1]) * (trans[1][0]) * (emitL[die.indexOf(input[i])]); //LL transition --> Loaded die to loaded die
      output[1][i] = (output[5][i-1]) * (trans[2][0]) * (emitL[die.indexOf(input[i])]); // FL transition
      output[2][i] = output[0][i];
      if(output[1][i] > output[0][i]) {
        path[0][i] = 1;
        output[2][i] = output[1][i];
      }
      
    
      output[3][i] = (output[2][i-1]) * (trans[1][1]) * (emitF[die.indexOf(input[i])]); //LF transition
      output[4][i] = (output[5][i-1]) * (trans[2][1]) * (emitF[die.indexOf(input[i])]); // FF transition
      output[5][i] = output[3][i];
      if(output[4][i] > output[3][i]) {
        path[1][i] = 1;
        output[5][i] = output[4][i];
      }
    }
    
  print2Array(output, 6, input.length);


    // traceback
    // To Get the viterbi path 
    int i = input.length -1;
    String vPath;
    int c;
    if(output[5][i] > output[2][i]){
      vPath = "1";
      c = path[1][i]; 
    } else {
      vPath = "0";
      c = path[0][i];
    } 
    for( i=input.length-2; i>=0; i--) {
     if(c == 0) {
        vPath = "0" + vPath;
      } else {
        vPath = "1" + vPath;
      }
      c = path[c][i];
    }
    System.out.println(vPath);

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
	    File f = new File("NC_000909.fna");
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
