// --------------------------------------------------------------
//  Evolve an HP mechanism to recover pyloric function best from 
//  a set of prescribed starting conditions
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

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
const int POPSIZE = 50;
const int GENS = 100;
const int trials = 10;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;
const double scaling_factor = 25; // boost to add to solutions that are fully pyloric

// Parameter variability modality only
//const int Repetitions = 10; 
//const int AnalysisReps = 100;

// Nervous system params
// const int N = 3;

// Plasticity parameters
const double SWR = 10;		// Max Window Size of Plastic Rule (in seconds now)
const double LBMIN = 0;
const double UBMIN = 0; 		// Plasticity Boundaries
const double LBMAX = 1;
const double UBMAX = 1; 		// Plasticity Boundaries 
const double BTMIN = 100.0;		// Bias Time Constant
const double BTMAX = 200.0;		// Bias Time Constant
// const double WTMIN = 40.0;		// Weight Time Constant
// const double WTMAX = 40.0;		// Weight Time Constant

int	VectSize = 8;

int num = 2; // Number of parameters being changed--temporary
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
	// Upper Bounds -- must be greater than lower bound
	for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), UBMIN, UBMAX);
		if (phen(k) < phen(k-2)){phen(k)=phen(k-2);}
		k++;
	}
    // Sliding Window -- changed to be time-based (gets rounded to the nearest stepsize in the SetSlidingWindow function)
    for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), 1, SWR);
		k++;
	}
}
ofstream pyloricoutput;

// ------------------------------------
// Recovery Fitness Function
// ------------------------------------
double HPFitnessFunction(TVector<double> &genotype, RandomState &rs){
    
    // Map genootype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
	
	// Create the agent
	CTRNN Agent(3);

	// Instantiate the nervous system
	char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Agent; 
	ifs.close();

	// Instantiate the HP mechanism
	Agent.SetHPPhenotype(phenotype,StepSize);

	double fitness = HPPerformance(Agent, scaling_factor);

    return fitness; //fitness averaged across all times it is taken
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

	BestIndividualsFile << trial << endl;
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
	BestIndividualsFile.open("bestinds.dat");
	for (int i=1;i<=trials;i++){
		long IDUM=-time(0);
		TSearch s(VectSize);

		// Configure the search
		s.SetRandomSeed(IDUM);
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

		s.SetEvaluationFunction(HPFitnessFunction);
		s.ExecuteSearch();

	// ifstream genefile("bestinds.dat");
	// int gen;
	// genefile >> gen;
	// TVector<double> genotype(1, VectSize);
	// genefile >> genotype;
	// cout << genotype << endl;
	// RandomState rs;
	// cout << HPFitnessFunction(genotype,rs) << endl;
		
	}
	Evolfile.close();
	BestIndividualsFile.close();
	pyloricoutput.close();

  return 0;
}
