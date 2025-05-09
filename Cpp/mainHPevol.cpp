// --------------------------------------------------------------
//  Evolve an HP mechanism to recover pyloric function best from 
//  a set of prescribed starting conditions
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"
#include "VectorMatrix.h"
#include <stdio.h>
#include <cstring>
#include <sys/stat.h> 

//#define PRINTOFILE

// Task params (using defined in the pyloric file)
const double TransientDuration = 250; //seconds without HP
const double PlasticDuration1 = 5000; //seconds allowing HP to act
const double PlasticDuration2 = 5000; //number of seconds to wait before testing again, to make sure not relying on precise timing
// const bool HPtest = true;       //does HP remain on during test (shouldn't matter if plasticity time constants are slow enough)
// const double StepSize = 0.25;
// const int TestSteps = TestDuration/StepSize; // in steps

// Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .1; //for detecting double periodicity

// EA params
const int POPSIZE = 25;
const int GENS = 100;
const int trials = 1;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;
const double pyloric_boost = 25; // boost to add to solutions that are fully pyloric
const bool seed_center_crossing = false;

const int num_optimization_genomes = 20; //number of circuits to test a generalist mechanism on (number listed in the file)

// Parameter variability modality only
//const int Repetitions = 10; 
//const int AnalysisReps = 100;

// Nervous system params
const int N = 3;
const bool shiftedrho_tf = true;
	//file from which to pull circuit genome
	// INDIVIDUAL IN EACH FOLDER MODE
const char circuitfname[] = "../pyloriccircuit.ns";
	// ONE INDIVIDUAL MODE
// const char circuitfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/92/pyloriccircuit.ns";
	// LIST OF INDIVIDUALS FOR THE GENERALIST MODE
// const char circuitfname[] = "../../../Pyloric CTRNN Genomes/optimizationsetengineered.dat";

// HP params
const int num = 2; //number of parameters you want to evolve to change under HP control
const int VectSize = num*4; //each parameter gets a time constant, an lower bound, a range, and a sliding window

// Plasticity parameter ranges
const double SWR = 10;		// Max Window Size of Plastic Rule (in seconds now)
const double LBMIN = 0;
const double UBMIN = 0; 		// OR range for range encoding
const double LBMAX = 1;
const double UBMAX = 1; 		// OR range for range encoding
const double BTMIN = 100.0;		// parameter Time Constant
const double BTMAX = 200.0;		// parameter Time Constant
// const double WTMIN = 40.0;		// Weight Time Constant
// const double WTMAX = 40.0;		// Weight Time Constant

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Bias Time-constants
	for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), BTMIN, BTMAX);
		k++;
	}
	// Lower Bounds
	for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), LBMIN, LBMAX);
		k++;
	}
	// Upper Bounds OR Range
	for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), UBMIN, UBMAX);
		// if (phen(k) < phen(k-num)){phen(k)=phen(k-num);} //clipping now happens in the Set functions
		k++;
	}
    // Sliding Window -- changed to be time-based (gets rounded to the nearest stepsize in the SetSlidingWindow function)
    for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), 0, SWR);
		k++;
	}
}

double HPPerformance(CTRNN &Agent, TMatrix<double> &ptlist, double pyloric_boost){
	double fitness = 0;
	int num_points = ptlist.RowUpperBound();
	// cout << "plane" << Agent.biases << endl;

	for (int i = 1; i <= num_points; i ++){
		for (int b=1;b<=num;b++){
			Agent.SetArbDParam(b,ptlist(i,b)); 
		}
		Agent.WindowReset();
		// cout << "init " << Agent.biases << endl;
		// cout << "parameters " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
		// Initialize the outputs at 0.5 for all neurons in the circuit
		Agent.RandomizeCircuitOutput(.5,.5);

		// Run the circuit for an initial transient; HP is off and fitness is not evaluated
		for (double t = StepSize; t <= TransientDuration; t += StepSize) {
			Agent.EulerStep(StepSize,false);
		}

		// Run the circuit for a period of time with HP so the paramters can change
		for (double t = StepSize; t<= PlasticDuration1; t+= StepSize){
			Agent.EulerStep(StepSize,true);  //set to only adapt biases, not weights
		}

		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";
		double fit = PyloricPerformance(Agent);
		// cout << "fitness" << fit;

		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+pyloric_boost;}

		// cout << "scaledfit" << fit;

		// cout << fit << endl;
		fitness += fit;

		// repeat process after another duration of time to ensure solution is stable
		for (double t = StepSize; t<= PlasticDuration2; t+= StepSize){
			Agent.EulerStep(StepSize,true);
		}

		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";
		fit = PyloricPerformance(Agent);

		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+pyloric_boost;}

		// cout << "fitness " << i << " " << fit << endl;
		fitness += fit;
		// cout << "final " << Agent.biases << endl;
		// cout << fit << endl << endl;
	}
	// cout << Agent.biases << endl;
    return fitness/(num_points*2);
}

// ------------------------------------
// Recovery Fitness Function
// ------------------------------------
//question for future lindsay: is this normalized by the fitness? probably not, so should we in order to compare solutions to different circuits?
double HPFitnessFunction(TVector<double> &genotype, TMatrix<double> &ptlist, RandomState &rs){
    // Map genotype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
	// cout << genotype << endl;
	// cout << "mapped";
	
	// Create the agent
	CTRNN Agent(3);

	// Instantiate the nervous system
    ifstream ifs;
    ifs.open(circuitfname);
    if (!ifs) {
        cerr << "File not found: " << circuitfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Agent; 
	ifs.close();
	Agent.ShiftedRho(shiftedrho_tf);
	// cout << "agent created" << endl;

	// Instantiate the HP mechanism
	
	// cout << Agent.PlasticityLB(1) << " " << Agent.PlasticityLB(2) << " " << Agent.PlasticityLB(3) << endl;
	// cout < < phenotype << endl;
	Agent.SetHPPhenotype(phenotype,StepSize,true); //range encoding active
	// ifstream HPf;
	// HPf.open("./Convenient HP Mechanisms/nullHP.dat");
	// Agent.SetHPPhenotype(HPf,StepSize,true);

	// cout << "plasticity set " << Agent.PlasticityLB(1) << " "  << Agent.PlasticityUB(1) << endl;
	double fitness = HPPerformance(Agent, ptlist, pyloric_boost);
	// cout << "fitness calculated" << endl;
    return fitness; //fitness averaged across all times it is taken
}

double HPGeneralistPerformance(CTRNN &Agent, ifstream &optimizationgenomesfile, TMatrix<double> &ptlist){
	//Here, Agent is already loaded with HP mechanism

	TVector<double> CTRNNphenotype;
	CTRNNphenotype.SetBounds(1,(2*N)+(N*N));

	TVector <double> performances(1,num_optimization_genomes);
	TVector <double> pyloric_fitnesses(1,num_optimization_genomes);

	for(int i=1;i<=num_optimization_genomes;i++){
		optimizationgenomesfile >> CTRNNphenotype;
		// cout << CTRNNphenotype << endl;
		optimizationgenomesfile >> pyloric_fitnesses[i]; //load in the fitness at the peak, so we can scale HP's performance

		CTRNNphenotype >> Agent;
		Agent.RandomizeCircuitOutput(.5,.5);
		//test the HP mechanism in each of the subspaces given by the evolved pyloric solutions
		
		performances[i] = HPPerformance(Agent, ptlist,pyloric_boost);
		// cout << performances[i] << endl << endl;
	}

	double max_pyloric_fitness = pyloric_fitnesses.Max();
	// cout << endl << max_pyloric_fitness << endl;
	double fitness = 0;
	//scale each performance by the peak (best possible performance relative to other solutiouns) 
	//and take the average
	for(int i=1;i<=num_optimization_genomes;i++){
		fitness += performances[i]/(pyloric_fitnesses[i]/max_pyloric_fitness);
	}
	// ifs.close();
	return fitness/num_optimization_genomes;
}

double HPGeneralistFitnessFunction(TVector<double> &genotype, TMatrix<double> &ptlist, RandomState &rs){
	// Map genotype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
	// cout << "mapped";
	
	// Create the agent
	CTRNN Agent(3);
	Agent.ShiftedRho(shiftedrho_tf);
	// cout << Agent.plasticneurons << endl;
	// cout << Agent.adaptweights << endl;
	// cout << "checkpoint 1" << endl;

	//Set HP parameters
	// cout << phenotype << endl;
	Agent.SetHPPhenotype(phenotype,StepSize,true); //using a phenotype vector, with range encoding
	// cout << Agent.l_boundary << " " << Agent.u_boundary << endl;
	// cout << "checkpoint 2" << endl;

	// Load in list of genomes in the optimization set
	ifstream ifs;
    ifs.open(circuitfname);
    if (!ifs) {
        cerr << "File not found: " << circuitfname << endl;
        exit(EXIT_FAILURE);
    }
	// cout << "checkpoint 3" << endl;

	double fit = HPGeneralistPerformance(Agent,ifs,ptlist);
	return fit;
}



// ------------------------------------
// Display functions
// ------------------------------------
ofstream Evolfile;
ofstream BestIndividualsFile;

int trial = 1;
void ResultsDisplay(TSearch &s)
{
	TVector<double> bestVector;
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);

	// Save the genotype of the best individual
	bestVector = s.BestIndividual();
	GenPhenMapping(bestVector, phenotype);

	// Reproduce which pars the HP mechanism has access to
	//LOCAL MODE
	// char plasticparsfname[] = "./plasticpars.dat";
	//SUPERCOMPUTER MODE
	char plasticparsfname[] = "../../../../plasticpars.dat";
  	ifstream plasticparsfile;
  	TVector<int> plasticitypars(1,N+(N*N));
  	plasticparsfile.open(plasticparsfname);
	if (!plasticparsfile) {
        cerr << "File not found: " << plasticparsfname << endl;
        exit(EXIT_FAILURE);
    }
  	for (int i = 1; i <= plasticitypars.UpperBound(); i ++){
    	plasticparsfile >> plasticitypars[i];
  	}
	// cout << plasticitypars << endl;
	// BestIndividualsFile << trial << endl;
	BestIndividualsFile << plasticitypars << endl;
	BestIndividualsFile << bestVector << endl << phenotype << endl;
	BestIndividualsFile << s.BestPerformance() << endl << endl;

	// cout << trial << "finished" << endl;

	trial ++;
}

void EvolutionaryRunDisplay(TSearch &s)
{
	
	//cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
	Evolfile << s.Generation() << " " << s.BestPerformance() << " " << s.AvgPerformance() << " " << s.PerfVariance() << endl;

	TVector<double> bestVector;
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);

	// Save the genotype of the best individual
	bestVector = s.BestIndividual();
	GenPhenMapping(bestVector, phenotype);

	Evolfile << phenotype << endl;
}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// SUPERCOMPUTER MODE
	Evolfile.open("./evol.dat");
	BestIndividualsFile.open("./bestind.dat");
	//SPECIFIC LOCAL MODE
	// char dirname[] = "./Specifically Evolved HP mechanisms/Every Circuit/92/0";
	// int result = mkdir(dirname,0755);
	// Evolfile.open("./Specifically Evolved HP mechanisms/Every Circuit/92/0/evol.dat");
	// BestIndividualsFile.open("./Specifically Evolved HP mechanisms/Every Circuit/92/0/bestind.dat");

	for (int i=1;i<=trials;i++){
		// cout << "trial " << i << endl;
		long randomseed = static_cast<long>(time(NULL));
		// long IDUM=-time(0);
		// cout << "seed initialized" << endl;
		// No idea why this wouldn't work
		// if (argc == 1)
		// 	randomseed += atoi(argv[1]);
		// cout << "seed set" << endl;
		TSearch s(VectSize);
		// cout << "Tsearch set up" << endl;

		// Configure the search
		s.SetRandomSeed(randomseed);
		//cout << IDUM << endl;
		s.SetSearchResultsDisplayFunction(ResultsDisplay);
		s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
		s.SetSelectionMode(RANK_BASED);
		s.SetReproductionMode(GENETIC_ALGORITHM);
		s.SetPopulationSize(POPSIZE);
		s.SetMaxGenerations(GENS);
		s.SetCrossoverProbability(CROSSPROB);
		s.SetCrossoverMode(UNIFORM);
		s.SetMutationVariance(MUTVAR);
		s.SetMaxExpectedOffspring(EXPECTED);
		s.SetElitistFraction(ELITISM);
		s.SetSearchConstraint(1);
		s.SetReEvaluationFlag(0); //  Parameter Variability Modality Only

		// GRID MODE
		int resolution = 5;
		TVector<double> par_vals(1,resolution);
		par_vals[1] = -16;
		par_vals[2] = -8;
		par_vals[3] = 0;
		par_vals[4] = 8;
		par_vals[5] = 16;

		int num_pts = pow(resolution,num);
		TMatrix<double> ptlist(1,num_pts,1,num);
		PointGrid(ptlist,par_vals);

		// RANDOM MODE
		// int num_pts = 50;
		
		// TMatrix<double> ptlist(1,num_pts,1,num);
		// for (int row = 1; row <= num_pts; row++){
		// 	for (int col = 1; col <= num; col++){
		// 		ptlist(row,col) = UniformRandom(-16,16);
		// 	}
		// }

		s.SetInitialPtsforEval(ptlist);
		s.SetEvaluationFunction(HPFitnessFunction);
		// cout << "TSearch configured" << endl;
		s.ExecuteSearch(seed_center_crossing);

		
	}
	Evolfile.close();
	BestIndividualsFile.close();

	// ifstream HPfile;
	// // ifstream Circuitfile;
	// HPfile.open("./converted33.dat");
	// // Circuitfile.open("./TestBestonDifferentSolutions/Pete/Pete.ns");

	// CTRNN Agent(3);
	// // Circuitfile >> Agent;
	// Agent.SetHPPhenotype(HPfile,StepSize,true);
	// // cout << Agent.l_boundary << " " << Agent.u_boundary << endl;

	// double fit = HPGeneralistPerformance(Agent,ptlist);
	// cout << fit;

	// HPfile.close();
	// Circuitfile.close();

  return 0;
}
