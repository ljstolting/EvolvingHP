// --------------------------------------------------------------
//  Evolve many circuits that satisfy the best HP mechanism as 
//  much as possible. How often are they pyloric? If relatively
//  often, this further confirms that there is something about the
//  averages that the HP mechanism has pulled out which is specific
//  to pyloricness in CTRNNs. Otherwise, it is only the case if
//  circuits start relatively close to/positioned in the right 
//  orientation to the closest circuit which satisfies HP and 
//  happens to be pyloric...
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "pyloric.h"
#include "random.h"

//#define PRINTOFILE

// Task params
// const double TransientDuration = 500; //in seconds
const double RunDuration = 25; //in seconds
// const double StepSize = 0.01;
const int RunSteps = RunDuration/StepSize; // in steps

// EA params
const int POPSIZE = 100;
const int GENS = 100;
const int trials = 2;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;
// const double ROC_weighting = 2000; // How heavily (relative to HPsatisfaction = 1) should maximizing the average |rate of change| be weighted?

// Nervous system params
const int N = 3;
const double WR = 16.0; 
const double BR = 16.0; //(WR*N)/2; //<-for allowing center crossing
const double TMIN = .1; 
const double TMAX = 2; 

int	VectSize = N*N + 2*N;

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), -BR, BR);
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				phen(k) = MapSearchParameter(gen(k), -WR, WR);
				k++;
			}
	}
}

// ----------------------
//  Fitness Function
// ------------------------

double HPsatisfaction(TVector<double> &phenotype){
	// HP mechanism
	char HPfname[] = "./33/bestind.dat";
	ifstream HPfile;
	HPfile.open(HPfname);

	CTRNN Circuit(N);
	// cout << Circuit.size << endl;
	phenotype >> Circuit; //just overwrote this operator to allow transfer from a phenotype vector
	// cout << "phenotype mapped " << Circuit.biases << endl;
	// cout << "phenotype mapped " << Circuit.biases << endl;
	bool range_encoding = true;
	Circuit.SetHPPhenotypebestind(HPfile,StepSize,range_encoding);
	// cout << "HP mapped" << Circuit.windowsize << endl;
	Circuit.RandomizeCircuitState(0,0);
	// cout << "phenotype mapped " << Circuit.biases << endl;
	TVector<double> acc(1,2);
	acc.FillContents(0);

	for(double t = StepSize; t <= TransientDuration; t += StepSize){
		Circuit.EulerStep(StepSize,false,false);
	}
	// cout << Circuit.outputs << endl;
	// cout << "Transient complete" << endl;
	for (double t=StepSize;t<=RunDuration;t+=StepSize){
		Circuit.EulerStepAvgsnoHP(StepSize);
		acc(1) += StepSize * Circuit.RtausBiases[1] * Circuit.rhos[1];
		acc(2) += StepSize * Circuit.RtausBiases[3] * Circuit.rhos[3];
		// cout << "single step complete" << endl;
		// cout << abs(acc(1)) << endl;
	}
	// cout << acc << " " << abs(acc(1))+abs(acc(2)) << endl;
	return 	1/(abs(acc(1))+abs(acc(2))); //no motivation for it to be small necessarily, just balanced. Is this reasonable?
}

double RateofChange(TVector<double> &phenotype){
	//returns the sum of the averages of the absolute values of the rates of change of the states for each neuron
	CTRNN Circuit(N);
	// cout << Circuit.size << endl;
	phenotype >> Circuit; //just overwrote this operator to allow transfer from a phenotype vector
	// cout << "phenotype mapped " << Circuit.biases << endl;
	// cout << "phenotype mapped " << Circuit.biases << endl;
	bool range_encoding = true;
	// cout << "HP mapped" << Circuit.windowsize << endl;
	Circuit.RandomizeCircuitState(0,0);
	for(double t = StepSize; t <= TransientDuration; t += StepSize){
		Circuit.EulerStep(StepSize,false,false);
	}
	// cout << Circuit.outputs << endl;
	// cout << "Transient complete" << endl;
	TVector<double> temp(1,N);
	TVector<double> cumulative_ROC(1,N);
	int step_counter = 0;
	for (double t=StepSize;t<=RunDuration;t+=StepSize){
		temp = Circuit.states;
		Circuit.EulerStep(StepSize,false,false);
		for (int i = 1; i <= N; i ++){
			cumulative_ROC[i] += abs(temp[i] - Circuit.NeuronState(i));
		}
		step_counter++;
	}
double fit = 0;
for (int i = 1; i <= N; i ++){
	fit += cumulative_ROC[i]/step_counter;
} 

if (fit <= .01) {fit = 0;} //some cutoff to define oscillation

return fit;

}

double ComboFitness(TVector<double> &genotype, RandomState &rs){
	//returns the weighted average of the avg absolute value of the rate of change and the HPsatisfaction function
	TVector<double> phenotype(1,VectSize);
	GenPhenMapping(genotype,phenotype);
	// cout << phenotype << endl;

	// double fitness = (HPsatisfaction(phenotype) + (ROC_weighting*RateofChange(phenotype)))/(1+ROC_weighting);
	double fitness  = HPsatisfaction(phenotype) * RateofChange(phenotype); //must be oscillating to get any fitness from the HP satisfaction (could make more of a threshold)
	return fitness;
}

double ROCFitonly(TVector<double> &genotype, RandomState &rs){
	//returns the weighted average of the avg absolute value of the rate of change and the HPsatisfaction function
	TVector<double> phenotype(1,VectSize);
	GenPhenMapping(genotype,phenotype);
	// cout << phenotype << endl;

	double fitness = RateofChange(phenotype);
	return fitness;
}

double HPFitonly(TVector<double> &genotype, RandomState &rs){
	//returns the weighted average of the avg absolute value of the rate of change and the HPsatisfaction function
	TVector<double> phenotype(1,VectSize);
	GenPhenMapping(genotype,phenotype);
	// cout << phenotype << endl;

	double fitness = HPsatisfaction(phenotype);
	return fitness;
}

double Pyloricfit(TVector<double> &phenotype){
	CTRNN Circuit(N);
	// cout << Circuit.size << endl;
	phenotype >> Circuit; //just overwrote this operator to allow transfer from a phenotype vector
	// cout << "phenotype mapped " << Circuit.biases << endl;
	// cout << "phenotype mapped " << Circuit.biases << endl;
	
	Circuit.RandomizeCircuitState(0,0);

	return PyloricPerfwTransient(Circuit); //hopefully this also prints the start and end steps of each neuron to determine which ordering criteria is breeched
}


// ------------------------------------
// Display functions
// ------------------------------------
ofstream Evolfile;
ofstream BestIndividualsFile;
ofstream trajfile;
ofstream ord_criteria_file;

void EvolutionaryRunDisplay(TSearch &s)
{
	
	Evolfile << s.Generation() << " " << s.BestPerformance() << " " << s.AvgPerformance() << " " << s.PerfVariance() << endl;
	// Evolfile << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}
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
	BestIndividualsFile << s.BestPerformance() << " " << HPsatisfaction(phenotype) << " " << RateofChange(phenotype) << endl << endl;

	cout << trial << "finished" << endl;

	trial ++;
}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{

	Evolfile.open("evolCTRNNfastcombo.dat");
	BestIndividualsFile.open("bestindsCTRNNfastcombo.dat");

	// trajfile.open("evolvedcircuit.dat");

	ord_criteria_file.open("orderingcriteriafastcombo.dat");

	for (int i=1;i<=trials;i++){
		long IDUM=-time(0);
		TSearch s(VectSize);

		#ifdef PRINTOFILE
		ofstream file;
		file.open("evol.dat");
		cout.rdbuf(file.rdbuf());
		#endif

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

		s.SetEvaluationFunction(ROCFitonly);
		s.ExecuteSearch(true);

		// ifstream genefile("best.gen.dat");
		// TVector<double> genotype(1, VectSize);
		// genefile >> genotype;

		//Gather trajectory of best evolved circuit
		TVector<double> bestgenotype(1, VectSize);
		TVector<double> bestphenotype(1, VectSize);
		bestgenotype = s.BestIndividual();
		GenPhenMapping(bestgenotype,bestphenotype);
		CTRNN Circuit(N);
		bestphenotype >> Circuit;
		Circuit.RandomizeCircuitOutput(.5,.5);

		// for(double t = StepSize; t <= TransientDuration; t += StepSize){
		// 	Circuit.EulerStep(StepSize,false,false);
		// 	trajfile << Circuit.outputs << endl;
		// }

		// for (double t=StepSize;t<=RunDuration;t+=StepSize){
		// 	Circuit.EulerStep(StepSize,false,false);
		// 	// cout << "single step complete" << endl;
		// 	trajfile << Circuit.outputs << endl;
		// }
		// trajfile << endl;

		//record the ordering criteria of the best evolved circuit
		OrderingRecord(Circuit,ord_criteria_file);
	}
	Evolfile.close();
	BestIndividualsFile.close();
	// trajfile.close();
	ord_criteria_file.close();
  return 0;
}
