// --------------------------------------------------------------
//  Evolve an HP mechanism to recover pyloric function best from 
//  a set of prescribed starting conditions
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"
#include "VectorMatrix.h"

//#define PRINTOFILE

// Task params (using defined in the pyloric file)
// const double TransientDuration = 250; //seconds without HP
// const double PlasticDuration1 = 2500; //seconds allowing HP to act
// const double PlasticDuration2 = 2500; //number of seconds to wait before testing again, to make sure not relying on precise timing
// const double TestDuration = 250; //maximum number of seconds allowed to test pyloric performance -- can be with HP still on
// const bool HPtest = true;       //does HP remain on during test (shouldn't matter if platicity time constants are slow enough)
// const double StepSize = 0.25;
// const int TestSteps = TestDuration/StepSize; // in steps

// Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .1; //for detecting double periodicity

// EA params
const int POPSIZE = 20;
const int GENS = 100;
const int trials = 1;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;
const double scaling_factor = 25; // boost to add to solutions that are fully pyloric

const int num_optimization_genomes = 10; //number of circuits to test a generalist mechanism on (number listed in the file)

// Parameter variability modality only
//const int Repetitions = 10; 
//const int AnalysisReps = 100;

// Nervous system params
const int N = 3;

// Plasticity parameter ranges
const double SWR = 10;		// Max Window Size of Plastic Rule (in seconds now)
const double LBMIN = 0;
const double UBMIN = 0; 		// Plasticity Boundaries
const double LBMAX = 1;
const double UBMAX = 1; 		// Plasticity Boundaries 
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
	for (int i = 1; i <= neuronschanging; i++) {
		phen(k) = MapSearchParameter(gen(k), LBMIN, LBMAX);
		k++;
	}
	// Upper Bounds OR Range
	for (int i = 1; i <= neuronschanging; i++) {
		phen(k) = MapSearchParameter(gen(k), UBMIN, UBMAX);
		// if (phen(k) < phen(k-num)){phen(k)=phen(k-num);} //clipping now happens in the Set functions
		k++;
	}
    // Sliding Window -- changed to be time-based (gets rounded to the nearest stepsize in the SetSlidingWindow function)
    for (int i = 1; i <= neuronschanging; i++) {
		phen(k) = MapSearchParameter(gen(k), 0, SWR);
		k++;
	}
}

double HPPerformance(CTRNN &Agent, TMatrix<double> &ptlist, double scaling_factor){
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
		for (int n=1;n<=Agent.CircuitSize();n++){Agent.SetNeuronState(n,0);}

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
		if (fit >= .3){fit = fit+scaling_factor;}

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
		if (fit >= .3){fit = fit+scaling_factor;}

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
double HPFitnessFunction(TVector<double> &genotype, TMatrix<double> &ptlist, RandomState &rs){
    // Map genotype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
	// cout << "mapped";
	
	// Create the agent
	CTRNN Agent(3);
	// cout << Agent.adaptbiases << endl;

	// Instantiate the nervous system
	char fname[] = "./Functioning2D/TestBestonDifferentSolutions/Pete/Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Agent; 
	ifs.close();

	// Instantiate the HP mechanism
	
	// cout << Agent.PlasticityLB(1) << " " << Agent.PlasticityLB(2) << " " << Agent.PlasticityLB(3) << endl;
	Agent.SetHPPhenotype(phenotype,StepSize,true); //range encoding active
	// cout << Agent.PlasticityLB(1) << " " << Agent.PlasticityLB(2) << " " << Agent.PlasticityLB(3) << endl;
	double fitness = HPPerformance(Agent, ptlist, scaling_factor);
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
		
		performances[i] = HPPerformance(Agent, ptlist,scaling_factor);
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
	// cout << Agent.plasticneurons << endl;
	// cout << Agent.adaptweights << endl;
	// cout << "checkpoint 1" << endl;

	//Set HP parameters
	// cout << phenotype << endl;
	Agent.SetHPPhenotype(phenotype,StepSize,true); //using a phenotype vector, with range encoding
	// cout << Agent.l_boundary << " " << Agent.u_boundary << endl;
	// cout << "checkpoint 2" << endl;

	// Load in list of genomes in the optimization set
	char fname[] = "../../../Pyloric CTRNN Genomes/optimizationset3.dat";
	ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
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
	char plasticparsfname[] = "../../plasticpars.dat";
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
	cout << plasticitypars << endl;
	// BestIndividualsFile << trial << endl;
	BestIndividualsFile << plasticitypars << endl;
	BestIndividualsFile << bestVector << endl << phenotype << endl;
	BestIndividualsFile << s.BestPerformance() << endl << endl;

	cout << trial << "finished" << endl;

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
	// Evolution condition
	Evolfile.open("evol.dat");
	BestIndividualsFile.open("bestind.dat");
	for (int i=1;i<=trials;i++){
		long randomseed = static_cast<long>(time(NULL));
		if (argc == 2)
			randomseed += atoi(argv[1]);
		// long IDUM=-time(0);
		TSearch s(VectSize);

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
		int resolution = 3;
		TVector<double> par_vals(1,resolution);
		par_vals[1] = -10;
		par_vals[2] = 0;
		par_vals[3] = 10;

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
		s.SetEvaluationFunction(HPGeneralistFitnessFunction);
		s.ExecuteSearch(false);

		
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
