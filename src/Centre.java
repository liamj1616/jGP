// MAKE SURE U UPDATE R H AND N TO MATCH THE DATASET!!!
// Change Dataset on line 133
// Change H on line 103
// Change R on line 105
// Change N on line 12

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Centre
{

	static int SEED;

	static int NUMBER_RUNS;
	static int POPULATION_SIZE;
	static int NUMBER_ISLANDS;
	static int NUMBER_MIGRATIONS;
	static int GENERATIONS;
	static int NUMBER_GENES;
	static double CROSSOVER_RATE;
	static double MUTATION_RATE;
	static int MUTATION_NUMBER;
	static int RING_DISTANCE_RATE;

	static int NUMBER_TRAINERS;
	static int NUMBER_PREDICTORS;
	static double PREDICTOR_SIZE;
	//static double PREDICTOR_GENS;

	static int SEGMENT;

	static Chromosome[][] curPop;
	static Chromosome[][] newPop;
	static Chromosome[] wholePop;

	static Thread[] islands;

	static Chromosome[] trainers;
	static double[][][] predictors;
	static double[][][] newPredictors;

	static double[][] allData;				// I added this for options pricing (this is the whole set of data, and theData is only 2 cols (for 1D data))
	static double[][] cashflowMatrix;		// I added this for options pricing (needed for knowing what was beign regressed and stuff)
	static boolean[][] didRegressionIndxs;	// I added this for options pricing (Needed it for knowing which rows were actually regressed
	static double[][] theData;

	static double[][] stats;

	static Random RNG;


	static double STRIKE_PRICE;
	static double R;
	static double H;
	static double DISCOUNT_RATE;

	static Chromosome[] eachRegressionTopFunctions;


	//static int[] testingPredictors = new int[122];

	public static void main(String[] args) throws InterruptedException
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		SEED = 1;
		NUMBER_RUNS = 1;
		POPULATION_SIZE = 101;		//per island, i.e. POPULATION_SIZE * NUMBER_ISLANDS		ALSO, needs to be odd for elite (I was lazy)
		NUMBER_ISLANDS = 7;
		NUMBER_MIGRATIONS = 10;//	100 seemed way to high --- overfit.
		GENERATIONS = 101000;//2500			//25000 works well
		NUMBER_GENES = 32;	//40 worked very well // = 16; // = 32;
		CROSSOVER_RATE = 80;		//%			//DONE THIS WAY BECAUSE NEXT INT IS APPARENTLY FASTER THAN DOUBLE (with rand)
		MUTATION_RATE = 10;			//%
		MUTATION_NUMBER = 2;
		RING_DISTANCE_RATE = 5;

		NUMBER_TRAINERS = 15;
		NUMBER_PREDICTORS = 20;
		///////////////////////////////
		PREDICTOR_SIZE = 0.1;		//rate based on theData size (NORMALLY 0.1, FOR COMBINATION, 0.01!!)
		//
		//PREDICTOR_GENS = 0.05;      //rate based on GENERATIONS (NOT GENES)

		/*
		 * args[0] = file name
		 * args[1] = which run are we doing
		 * args[2] = how many tries per regression are we doing
		 * args[3] = STRIKE_PRICE
		 * args[4] = H
		 * args[5] = R
		 *
		 */

		/// stuff for options
		STRIKE_PRICE = 40.00;		// 55.00
		// H = Double.parseDouble(args[4]);				// 0.0027397260273972603
		H = 0.0833333;
		// R = Double.parseDouble(args[5]);				// 0.02
		R = 0.0601;
		DISCOUNT_RATE = Math.exp(-1*R*H);
		/// stuff for options



		///////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//String fileName = "/scratch/m/mdaley/jhughe54/outsOPT/";
		String fileName = "./outsOPT/";
		for (int i = 0; i < args.length; i++)
		{
			if(i != 0) fileName = fileName + "_";
			fileName = fileName + args[i];
		}


		new File(fileName).mkdir();
		fileName = fileName + "/";

		stats = new double[NUMBER_RUNS][];
		//tried GAMB, MOTOR, WM
		//theData = Reader.readData("../../../../WCCI2015/Walking/DATA/", "CH-FR-LACC-RAW-RM-5hz-1000-15000.csv");
		//theData = Reader.readData("../../optData/", "M1000_N12_S10.0_h0.01_sigma0.2_r0.05.csv");
		//theData = Reader.readData("../roiDATA/", "MOTOR_100307_2_L21_Z.csv");		//CHANGE THE PREDICTOR SIZE IF USING THIS!!!
		//theData = Reader.readData("/home/m/mdaley/jhughe54/optData/", args[0]+".csv");
		allData = Reader.readData("./data/", "generatedData30.csv");

		//Language.setVarSymbolsNumbers(theData[0].length - 1);
		//Language.setVarSymbols(Reader.readVarNames("./", "variableNames2.csv"));
		String[] varLabels = {"x"};
		Language.setVarSymbols(varLabels);



		cashflowMatrix = startCFM();
		didRegressionIndxs = new boolean[allData.length][allData[0].length];
		eachRegressionTopFunctions = new Chromosome[allData[0].length -1 - 1];		// Remember, fence posts and ignore (0-1), so if 12 time points, we only need 10 regressions (1-2, 2-3, 3-4, 4-5, etc.) *skip 0->1*


		/*
		for(int i = 0 ; i < allData.length ; i++)
		{
			System.out.println(allData[i][allData[0].length-1-1] + ", " + allData[i][allData[0].length-1] + "\t\t" + theData[i][0] + "," + theData[i][1]);
		}
		*/

		//theData = genDataForRegression(10);

		// loop for however many regressions we need to do
		int sValue = 11;
		double totalError = 0;

		for(int n = allData[0].length -1 -1; n >= 1 ; n--)	// counts down from n-1-1  to 1 (because if the data is 12, we would start at 11 and go to 1 (because we ignore time 0 apparently!?)
		{
			System.out.println("(" + (n) + " -> " + (n+1) +")");
			theData = genDataForRegression(n);		// change this to loops thing.
			if(theData.length == 0)
			{
				System.out.println("\t\t\tHELLO!");
				Chromosome dummyChromosome = new Chromosome(NUMBER_GENES, RNG);		// dummy chromosome when there is nothing to regress
				dummyChromosome.geneArray[NUMBER_GENES - 1] = new Gene('c', Double.MIN_VALUE);	// Just set the root to 0.0
				eachRegressionTopFunctions[n - 1] = dummyChromosome;
				continue;															// Just go to next
			}



			//TEST AREA

			ThreadKiller tKiller = new ThreadKiller();

			//for(int i = sValue ; i < sValue+NUMBER_RUNS ; i++)
			for(int i = 0 ; i < 0+NUMBER_RUNS ; i++)
			{
				System.out.println("\t" + i);
				//SEED = i;
				//SEED = 1;
				SEED = (int)System.currentTimeMillis();
				RNG = new Random(SEED);
				//int noLearnCount = 0;
				//double lastBestFitness = Double.MAX_VALUE;
				generatePopulations();

				Chromosome bestOverAllMigrations = new Chromosome(NUMBER_GENES, RNG);
				bestOverAllMigrations.fitnessValue = Double.MAX_VALUE;

				trainers = new Chromosome[NUMBER_TRAINERS];
				predictors = new double[NUMBER_PREDICTORS][][];
				newPredictors = new double[NUMBER_PREDICTORS][][];
				combinePops();
				TrainerWorker.initialize(trainers, wholePop, RNG);
				PredictorWorker.initialize(predictors, Math.max((int)(theData.length * PREDICTOR_SIZE), 1), theData, RNG);			// if # predictors are less than 0, then make it 1... This is hackish though. Think of something that makes more sense.
				PredictorWorker.initialize(newPredictors,  Math.max((int)(theData.length * PREDICTOR_SIZE), 1), theData, RNG);


				//loads the save state if the segment isn't zero
				//if(SEGMENT > 0)
				//{
				//StateSaveRestore.restoreStateSer("/scratch/m/mdaley/jhughe54/segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT - 1, bestOverAllMigrations, curPop, trainers, predictors);
				//	StateSaveRestore.restoreStateSer("./segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT - 1, bestOverAllMigrations, curPop, trainers, predictors);
				//}

				for(int j = 0 ; j < NUMBER_MIGRATIONS ; j++)
				{
					//System.out.println(j);
					tKiller.setDie(false);		//REVIVEEE
					Thread pool = new Thread(new PredictorPool(NUMBER_PREDICTORS, CROSSOVER_RATE, MUTATION_RATE, MUTATION_NUMBER, RNG, predictors, newPredictors, trainers, theData, tKiller));
					pool.start();
					for(int k = 0 ; k < islands.length ; k++)
					{
						islands[k] = new Thread(new Island(POPULATION_SIZE, GENERATIONS, NUMBER_MIGRATIONS, CROSSOVER_RATE, MUTATION_RATE, MUTATION_NUMBER, RING_DISTANCE_RATE, RNG, curPop[k], newPop[k], theData));
						islands[k].start();
					}
					for(int k = 0 ; k < islands.length ; k++)
					{
						islands[k].join();
					}
					tKiller.setDie(true);		//KILLER
					pool.join();	/////////////////////////////////////
					combinePops();

					Evaluator.calcAllFitnessMSE(wholePop, theData);
					Chromosome curBestResult = Evaluator.findBestChromosomeMin(wholePop);
					if(curBestResult.fitnessValue < bestOverAllMigrations.fitnessValue)
					{
						bestOverAllMigrations = curBestResult;
					}
					//System.out.println(bestOverAllMigrations.fitnessValue);
					/*
					for(int w = 0 ; w < wholePop.length ; w++)
					{
						System.out.print(wholePop[w].fitnessValue + ",");
					}
					System.out.println();
					System.out.println();
					*/

					TrainerWorker.newTrainers(trainers, wholePop, predictors);		////////////////////////////////
					shuffle(wholePop);
					separatePops();

					//System.out.println("\t\t\t\t " + PredictorWorker.evaluatePredictorFitness(PredictorWorker.bestPredictor, theData, trainers));
					//for(int k = 0 ; k < PredictorWorker.bestPredictor.length ; k++)
					//{
					//System.out.print("\t\t" + PredictorWorker.bestPredictor[k][0]);
					//testingPredictors[(int)((PredictorWorker.bestPredictor[k][0] + 3) * 20)]++;
					//}
					//System.out.println();

				}

				//
				//StateSaveRestore.saveStateSer("/scratch/m/mdaley/jhughe54/segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT, bestOverAllMigrations, curPop, trainers, predictors);
				//StateSaveRestore.saveStateSer("./segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT, bestOverAllMigrations, curPop, trainers, predictors);
				//StateSaveRestore.saveStateStr("./", SEGMENT, bestOverAllMigrations, curPop, trainers, predictors);
				//StateSaveRestore.deleteStateSer("/scratch/m/mdaley/jhughe54/segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT - 1);
				//StateSaveRestore.deleteStateSer("./segmentsOPT/" + args[0] + "_" + args[1] + "_", SEGMENT - 1);
				//

				combinePops();
				Evaluator.calcAllFitnessMSE(wholePop, theData);
				Chromosome curBestResult = Evaluator.findBestChromosomeMin(wholePop);
				if(curBestResult.fitnessValue < bestOverAllMigrations.fitnessValue)
				{
					bestOverAllMigrations = curBestResult;
				}
				//stats[i] = new double[3];
				//stats[i][0] = bestResult.fitnessValue;
				//stats[i][1] = Evaluator.calcAvgFitness(wholePop);
				//stats[i][2] = Evaluator.calcStandardDeviation(wholePop);


				// will need to figure out how to make printer work now
				//Printer.printLine(fileName + i, bestOverAllMigrations);
				//Printer.printMulti(fileName + i, bestOverAllMigrations);
				//Printer.printCurStat(fileName + i, bestOverAllMigrations);			//quick fix... I wish I did this a lot earlier (screw you automatic windows updates)
				System.out.println("\t\t" + bestOverAllMigrations.fitnessValue);
				totalError += bestOverAllMigrations.fitnessValue;


				//for(int k = 0 ; k < testingPredictors.length ; k++)
				//{
				//	System.out.println((((double)k/20)-3) + "\t" + testingPredictors[k]);
				//}
				if(eachRegressionTopFunctions[n - 1] == null || eachRegressionTopFunctions[n - 1].fitnessValue > bestOverAllMigrations.fitnessValue)
				{
					eachRegressionTopFunctions[n - 1] = new Chromosome(bestOverAllMigrations);
				}
			}
			//System.out.println(eachRegressionTopFunctions[n-1].fitnessValue);
			updateCFM(n, theData);
		}
		// Printer.printStats(fileName, stats);
		// need to update fileName here...
		Printer.printLineAllRegressors(fileName, sValue, eachRegressionTopFunctions);
		double totalValuation = 0.0;
		for (int i = 0; i < cashflowMatrix.length; i++){
			for (int j = 1; j < cashflowMatrix[i].length; j++){
				if (cashflowMatrix[i][j] != 0){
					totalValuation += cashflowMatrix[i][j] * Math.pow(DISCOUNT_RATE, j);
					break;
				}
			}
		}
		System.out.println(totalValuation / cashflowMatrix.length);
		System.out.println(totalError / sValue);
	}



	public static double[][] startCFM()
	{
		double[][] cfm = new double[allData.length][allData[0].length];	// we ignore time 0
		for(int i = 0 ; i < allData.length ; i++)
		{
			cfm[i][allData[0].length -1] = Math.max(0, STRIKE_PRICE - allData[i][allData[0].length-1]);		// I BELIEVE THIS BIT IS WHY I HAVE IT'S SET TO PUT (switch strike and data in minus)
		}
		return cfm;
	}

	public static double[][] genDataForRegression(int index)
	{
		double[][] regData;
		List<Double> tmpHolderX = new ArrayList<Double>();
		List<Double> tmpHolderY = new ArrayList<Double>();

		// This bit sets what should be in the X and Ys (including money discount)
		for(int i = 0 ; i < allData.length ; i++)
		{
			if(STRIKE_PRICE - allData[i][index] > 0)		// If it's ``in the money''
			{
				didRegressionIndxs[i][index] = true;
				tmpHolderX.add(allData[i][index]);

				// Go back towards end until we find the price that exists
				for(int j = index + 1 ; j < allData[0].length ; j++)
				{
					if(cashflowMatrix[i][j] > 0)
					{
						tmpHolderY.add(cashflowMatrix[i][j] * Math.pow(DISCOUNT_RATE, j - index));
						break;														// overkill --- If we ever find it, STOP
					}
					else
					{
						if(j == allData[0].length - 1)		// If we don't find anything, then it's zer0
						{
							tmpHolderY.add(0.00);
						}
					}
				}
			}
			else
			{
				didRegressionIndxs[i][index] = false;
			}
		}

		// This big just puts the data into a 2D array for theData
		regData = new double[tmpHolderX.size()][2];

		for(int i = 0 ; i < tmpHolderX.size() ; i++)
		{
			regData[i][0] = tmpHolderX.get(i);
			regData[i][1] = tmpHolderY.get(i);
		}

		return regData;
	}

	public static void updateCFM(int index, double[][] regData)
	{
		int count = 0;		// this count is for going through the
		double continuationPrice;
		double[] X = new double[1];     // has to be an array because that's what the evaluator function takes
		for(int i = 0 ; i < allData.length ; i++)
		{
			if(didRegressionIndxs[i][index])
			{
				X[0] = regData[count][0];
				continuationPrice = Evaluator.evaluateChromosome(eachRegressionTopFunctions[index - 1].numGenes - 1, X, eachRegressionTopFunctions[index - 1]);		// -1 because of of how we're regressing

				// WE COULD TOTALLY CALCULATE ERROR HERE!!!!!!
				// if we're IN THE MONEY!
				if (STRIKE_PRICE - X[0] >= continuationPrice)									// I BELIEVE THIS BIT IS WHY I HAVE IT'S SET TO PUT (switch strike and data in minus)
				{
					cashflowMatrix[i][index] = STRIKE_PRICE - X[0];								// update cfM
					// sets everything past should now be zer0
					for(int j = index + 1 ; j < cashflowMatrix[0].length ; j++)
					{
						cashflowMatrix[i][j] = 0.00;
					}

				}
				else		// if it's not, then make it 0, because we lose here
				{
					cashflowMatrix[i][index] = 0.00;
				}
				count++;

			}
			else		// if we didn't do anything, set CFM to 0.00... I really don't know if this is necessary.
			{
				cashflowMatrix[i][index] = 0.00;
			}
		}
		// we could return average error here if we wanted... just sayin'
	}

	// Could delete this later
	public static void makeSpecificChromosomesForTesting()
	{
		Random bleh = new Random(1234);
		eachRegressionTopFunctions[1] = new Chromosome(10, bleh);
		eachRegressionTopFunctions[1].geneArray[9] = new Gene('b', "+", 8, 0);
		eachRegressionTopFunctions[1].geneArray[8] = new Gene('b', "+", 7, 6);
		eachRegressionTopFunctions[1].geneArray[7] = new Gene('b', "*", 3, 1);
		eachRegressionTopFunctions[1].geneArray[6] = new Gene('b', "*", 5, 2);
		eachRegressionTopFunctions[1].geneArray[5] = new Gene('b', "*", 4, 3);
		eachRegressionTopFunctions[1].geneArray[4] = new Gene('v', 0);
		eachRegressionTopFunctions[1].geneArray[3] = new Gene('v', 0);		// I think I didn't need 2 Xs. Whatever
		eachRegressionTopFunctions[1].geneArray[2] = new Gene('c', -1.813);
		eachRegressionTopFunctions[1].geneArray[1] = new Gene('c', 2.983);
		eachRegressionTopFunctions[1].geneArray[0] = new Gene('c', -1.070);

		eachRegressionTopFunctions[0] = new Chromosome(10,bleh);
		eachRegressionTopFunctions[0].geneArray[9] = new Gene('b', "+", 8, 0);
		eachRegressionTopFunctions[0].geneArray[8] = new Gene('b', "+", 7, 6);
		eachRegressionTopFunctions[0].geneArray[7] = new Gene('b', "*", 3, 1);
		eachRegressionTopFunctions[0].geneArray[6] = new Gene('b', "*", 5, 2);
		eachRegressionTopFunctions[0].geneArray[5] = new Gene('b', "*", 4, 3);
		eachRegressionTopFunctions[0].geneArray[4] = new Gene('v', 0);
		eachRegressionTopFunctions[0].geneArray[3] = new Gene('v', 0);		// I think I didn't need 2 Xs. Whatever :p
		eachRegressionTopFunctions[0].geneArray[2] = new Gene('c', 1.356);
		eachRegressionTopFunctions[0].geneArray[1] = new Gene('c', -3.335);
		eachRegressionTopFunctions[0].geneArray[0] = new Gene('c', 2.038);
	}

	public static void generatePopulations()
	{
		curPop = new Chromosome[NUMBER_ISLANDS][];
		newPop = new Chromosome[NUMBER_ISLANDS][];

		wholePop = new Chromosome[POPULATION_SIZE * NUMBER_ISLANDS];

		islands = new Thread[NUMBER_ISLANDS];

		for(int i = 0 ; i < NUMBER_ISLANDS ; i++)
		{
			curPop[i] = new Chromosome[POPULATION_SIZE];
			newPop[i] = new Chromosome[POPULATION_SIZE];

			for(int j = 0 ; j < curPop[i].length ; j++)
			{
				curPop[i][j] = new Chromosome(NUMBER_GENES, RNG);
				curPop[i][j].genChromosome();
				newPop[i][j] = new Chromosome(NUMBER_GENES, RNG);
				newPop[i][j].genChromosome();
			}
		}
	}

	public static void combinePops()
	{
		int indexCounter = 0;

		for(int i = 0 ; i < NUMBER_ISLANDS ; i++)
		{
			for(int j = 0 ; j < curPop[i].length ; j++)
			{
				wholePop[indexCounter++] = curPop[i][j];
			}
		}
	}

	public static void separatePops()
	{
		int indexCounter = 0;

		for(int i = 0 ; i < NUMBER_ISLANDS ; i++)
		{
			for(int j = 0 ; j < curPop[i].length ; j++)
			{
				curPop[i][j] = wholePop[indexCounter++];
			}
		}
	}

	public static void shuffle(Chromosome[] inA)
	{
		Chromosome temp;
		int i = inA.length - 1;
		int j;

		while(i > 1)
		{
			j = RNG.nextInt(i+1);
			temp = inA[i];
			inA[i] = inA[j];
			inA[j] = temp;
			i--;
		}
	}





}


//
//import java.util.Random;
//
//public class Centre {
//
//	static int SEED;
//
//	static double STRIKE_PRICE;
//
//	static int NUMBER_RUNS;
//	static int POPULATION_SIZE;
//	static int NUMBER_ISLANDS;
//	static int NUMBER_MIGRATIONS;
//	static int GENERATIONS;
//	static int NUMBER_GENES;
//	static double CROSSOVER_RATE;
//	static double MUTATION_RATE;
//	static int MUTATION_NUMBER;
//	static int RING_DISTANCE_RATE;
//
//	static int NUMBER_TRAINERS;
//	static int NUMBER_PREDICTORS;
//	static double PREDICTOR_SIZE;
//	// static double PREDICTOR_GENS;
//
//	static int SEGMENT;
//
//	static Chromosome[][] curPop;
//	static Chromosome[][] newPop;
//	static Chromosome[] wholePop;
//
//	static Thread[] islands;
//
//	static Chromosome[] trainers;
//	static double[][][] predictors;
//	static double[][][] newPredictors;
//
//	static double[][] intermediaryData;
//	static double[][] theData;
//
//	static double[][] stats;
//
//	static Random RNG;
//
//	// static int[] testingPredictors = new int[122];
//
//	public static void main(String[] args) throws InterruptedException {
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//		STRIKE_PRICE = 40;
//
//		SEED = 1;
//		NUMBER_RUNS = 1;
//		POPULATION_SIZE = 101; // per island, i.e. POPULATION_SIZE * NUMBER_ISLANDS ALSO, needs to be odd for
//								// elite (I was lazy)
//		NUMBER_ISLANDS = 7;
//		NUMBER_MIGRATIONS = 100;// (CHANGE ME) //10000 was amazing (with 2500gens --- otherwise, wait like a
//								// day....)
//		GENERATIONS = 10100;//
//		NUMBER_GENES = 40; // 40 worked well // = 16; // = 32;
//		CROSSOVER_RATE = 80; // % //DONE THIS WAY BECAUSE NEXT INT IS APPARENTLY FASTER THAN DOUBLE (with
//								// rand)
//		MUTATION_RATE = 10; // %
//		MUTATION_NUMBER = 2;
//		RING_DISTANCE_RATE = 5;
//
//		NUMBER_TRAINERS = 8;
//		NUMBER_PREDICTORS = 10;
//		///////////////////////////////
//		PREDICTOR_SIZE = 0.25; // rate based on theData size (NORMALLY 0.1, FOR COMBINATION, 0.01!!)
//		//
//		// PREDICTOR_GENS = 0.05; //rate based on GENERATIONS (NOT GENES)
//
//		SEGMENT = 0;
//
//		///////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//		/*
//		 * String fileName = "outs/";
//		 * for (int i = 0; i < args.length; i++)
//		 * {
//		 * if(i != 0) fileName = fileName + "_";
//		 * fileName = fileName + args[i];
//		 * }
//		 *
//
//			new File(fileName).mkdir();
//			fileName = fileName + "/";
//		 */
//		 stats = new double[NUMBER_RUNS][];
//		 //tried GAMB, MOTOR, WM
//		 //theData = Reader.readData("../../../../WCCI2015/Walking/DATA/", "CH-FR-LACC-RAW-RM-5hz-1000-15000.csv");
//		 //theData = Reader.readData("../roiData/", "MOTOR_100307_3_L33.csv");
//		 // theData = Reader.readData("../roiDATA/", "MOTOR_100307_2_L21_Z.csv");
//		 //CHANGE THE PREDICTOR SIZE IF USING THIS!!!
//
//		String fileName = "./data/generatedData17.csv";
//		intermediaryData = Reader.readData("./data/", "generatedData17.csv");
//
//		theData = new double[intermediaryData.length][intermediaryData[0].length];
//		boolean flag = false;
//		int theDataLength = 0;
//		for (int i = 0; i < intermediaryData.length; i++){
//			flag = false;
//			for (int j = 1; j < intermediaryData[i].length; j++){
//				if (intermediaryData[i][j] < STRIKE_PRICE){
//					flag = true;
//					break;
//				}
//			}
//			if (flag == true){
//				theData[theDataLength] = intermediaryData[i];
//				theDataLength++;
//			}
//		}
//
//		System.out.println(theData.length);
//		System.out.println(theData[0].length);
//
//
//		// Sets the number of variables
//		Language.setVarSymbolsNumbers(theData[0].length - 1);
//
//
//		// TEST AREA
//		RNG = new Random(SEED);
//
//		ThreadKiller tKiller = new ThreadKiller();
//		for (int i = 0; i < NUMBER_RUNS; i++) {
//			System.out.println("Run number " + i);
//			SEED = (int) System.currentTimeMillis();
//		 	RNG = new Random(SEED);
//			 // int noLearnCount = 0;
//			// double lastBestFitness = Double.MAX_VALUE;
//		 	generatePopulations();
//
//			Chromosome bestOverAllMigrations = new Chromosome(NUMBER_GENES, RNG);
//			bestOverAllMigrations.fitnessValue = Double.MAX_VALUE;
//
//			trainers = new Chromosome[NUMBER_TRAINERS];
//			predictors = new double[NUMBER_PREDICTORS][][];
//			newPredictors = new double[NUMBER_PREDICTORS][][];
//			combinePops();
//
//			TrainerWorker.initialize(trainers, wholePop, RNG);
//			PredictorWorker.initialize(predictors, (int) (theData.length * PREDICTOR_SIZE), theData, RNG);
//			PredictorWorker.initialize(newPredictors, (int) (theData.length * PREDICTOR_SIZE), theData, RNG);
//
//			// loads the save state if the segment isn't zero
//			if (SEGMENT > 0) {
//				StateSaveRestore.restoreStateSer(
//						"./",
//						SEGMENT - 1,
//						bestOverAllMigrations,
//						curPop, trainers, predictors
//				);
//			}
//
//			for (int j = 0; j < NUMBER_MIGRATIONS; j++) {
//				tKiller.setDie(false); // REVIVEEE
//				Thread pool = new Thread(new PredictorPool(
//						NUMBER_PREDICTORS,
//						CROSSOVER_RATE,
//						MUTATION_RATE,
//						MUTATION_NUMBER,
//						RNG,
//						predictors,
//						newPredictors,
//						trainers,
//						theData,
//						tKiller
//				));
//
//				pool.start();
//
//				for (int k = 0; k < islands.length; k++) {
//					islands[k] = new Thread(new Island(
//							POPULATION_SIZE,
//							GENERATIONS,
//							NUMBER_MIGRATIONS,
//							CROSSOVER_RATE,
//							MUTATION_RATE,
//							MUTATION_NUMBER,
//							RING_DISTANCE_RATE,
//							RNG,
//							curPop[k],
//							newPop[k],
//							theData
//					));
//					islands[k].start();
//				}
//
//				for (int k = 0; k < islands.length; k++) {
//					islands[k].join();
//				}
//
//				tKiller.setDie(true); // KILLER
//				pool.join(); /////////////////////////////////////
//				combinePops();
//
//				Evaluator.calcAllFitnessMSE(wholePop, theData);
//				Chromosome curBestResult = Evaluator.findBestChromosomeMin(wholePop);
//
//				if (curBestResult.fitnessValue < bestOverAllMigrations.fitnessValue) {
//					bestOverAllMigrations = curBestResult;
//				}
//				// System.out.println(bestOverAllMigrations.fitnessValue);
//				/*
//				 * for(int w = 0 ; w < wholePop.length ; w++)
//				 * {
//				 * System.out.print(wholePop[w].fitnessValue + ",");
//				 * }
//				 * System.out.println();
//				 * System.out.println();
//				 */
//
//				// DELETE THIS TOO!!!!
//				TrainerWorker.newTrainers(trainers, wholePop, predictors);
//				////////////////////////////////
//				shuffle(wholePop);
//				separatePops();
//
//				// System.out.println("\t\t\t\t " +
//				// PredictorWorker.evaluatePredictorFitness(PredictorWorker.bestPredictor,
//				// theData, trainers));
//				// for(int k = 0 ; k < PredictorWorker.bestPredictor.length ; k++)
//				// {
//				// System.out.print("\t\t" + PredictorWorker.bestPredictor[k][0]);
//				// testingPredictors[(int)((PredictorWorker.bestPredictor[k][0] + 3) * 20)]++;
//				// }
//				// System.out.println();
//			}
//
//			StateSaveRestore.saveStateSer("./", SEGMENT, bestOverAllMigrations, curPop, trainers, predictors);
//			// StateSaveRestore.saveStateStr("./", SEGMENT, bestOverAllMigrations,curPop, trainers, predictors);
//			StateSaveRestore.deleteStateSer("./", SEGMENT - 1);
//
//
//			// stats[i] = new double[3];
//			// stats[i][0] = bestResult.fitnessValue;
//			// stats[i][1] = Evaluator.calcAvgFitness(wholePop);
//			// stats[i][2] = Evaluator.calcStandardDeviation(wholePop);
//
//			Printer.printLine(fileName + i, bestOverAllMigrations);
//			Printer.printMulti(fileName + i, bestOverAllMigrations);
//			Printer.printCurStat(fileName + i, bestOverAllMigrations); // quick fix... I wish I did this a lot earlier (screw you automatic windows updates)
//			System.out.println("\t\t" + bestOverAllMigrations.fitnessValue);
//
//			// for(int k = 0 ; k < testingPredictors.length ; k++)
//			// {
//			// System.out.println((((double)k/20)-3) + "\t" + testingPredictors[k]);
//			// }
//		}
//
//		// Printer.printStats(fileName, stats);
//
//	}
//
//	public static void generatePopulations() {
//		curPop = new Chromosome[NUMBER_ISLANDS][];
//		newPop = new Chromosome[NUMBER_ISLANDS][];
//
//		wholePop = new Chromosome[POPULATION_SIZE * NUMBER_ISLANDS];
//
//		islands = new Thread[NUMBER_ISLANDS];
//
//		for (int i = 0; i < NUMBER_ISLANDS; i++) {
//			curPop[i] = new Chromosome[POPULATION_SIZE];
//			newPop[i] = new Chromosome[POPULATION_SIZE];
//
//			for (int j = 0; j < curPop[i].length; j++) {
//				curPop[i][j] = new Chromosome(NUMBER_GENES, RNG);
//				curPop[i][j].genChromosome();
//				newPop[i][j] = new Chromosome(NUMBER_GENES, RNG);
//				newPop[i][j].genChromosome();
//			}
//		}
//	}
//
//	public static void combinePops() {
//		int indexCounter = 0;
//
//		for (int i = 0; i < NUMBER_ISLANDS; i++) {
//			for (int j = 0; j < curPop[i].length; j++) {
//				wholePop[indexCounter++] = curPop[i][j];
//			}
//		}
//	}
//
//	public static void separatePops() {
//		int indexCounter = 0;
//
//		for (int i = 0; i < NUMBER_ISLANDS; i++) {
//			for (int j = 0; j < curPop[i].length; j++) {
//				curPop[i][j] = wholePop[indexCounter++];
//			}
//		}
//	}
//
//	public static void shuffle(Chromosome[] inA) {
//		Chromosome temp;
//		int i = inA.length - 1;
//		int j;
//
//		while (i > 1) {
//			j = RNG.nextInt(i + 1);
//			temp = inA[i];
//			inA[i] = inA[j];
//			inA[j] = temp;
//			i--;
//		}
//	}
//
//}
