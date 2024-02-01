// --------------------------------------------------------------
//  Evolve an HP mechanism to recover pyloric function best from 
//  a set of prescribed starting conditions
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 500; //seconds without HP
const double PlasticDuration = 2000; //seconds allowing HP to act
const double TestDuration = 250; //maximum number of seconds allowed to test pyloric performance -- can be with HP still on
const bool HPtest = true;       //does HP remain on during test (shouldn't matter if platicity time constants are slow enough)
const double StepSize = 0.01;
const int TestSteps = TestDuration/StepSize; // in steps

// Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .05; //for detecting double periodicity

// EA params
const int POPSIZE = 50;
const int GENS = 250;
const int trials = 10;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;

// Parameter variability modality only
//const int Repetitions = 10; 
//const int AnalysisReps = 100;

// Nervous system params
const int N = 3;

// Plasticity parameters
const double SWR = 120;		// Max Window Size of Plastic Rule (in steps size) (so 1 is no window)
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
    // Sliding Window -- must be an integer, so round to nearest integer
    for (int i = 1; i <= num; i++) {
		phen(k) = MapSearchParameter(gen(k), 1, SWR);
        phen(k) = round(phen(k));
		k++;
	}
}

// ------------------------------------
// Pyloric-like performance
// ------------------------------------
double PyloricPerformance(CTRNN Agent)
{
    TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);
	double fitness = 0.0;

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness -- HP is either left on or turned OFF depending on specification.

	// Run the circuit until you identify 3 PD starts (capping 2 full cycles), keeping track of whether each other neuron crossed the threshold or not
	int tstep = 0;
	double time = StepSize;
	int PDstartcount = 0;
	TVector<int> PDstarts(1,3);
	PDstarts.FillContents(0);
	while (time <= TestDuration && PDstartcount < 3) {
		tstep += 1;
		for (int i = 1; i <= N; i += 1) {
			OutputHistory[tstep][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,HPtest,HPtest);

		//Check for PD start
		if (OutputHistory[tstep][3] < burstthreshold && Agent.NeuronOutput(3) > burstthreshold){
			PDstarts[3-PDstartcount] = tstep;
			PDstartcount += 1;
			//cout << "PDstarts";
		}
	}

	int criteriamet = 0;
	for (int i = 1; i <= N; i += 1) {
		// SHORT HAND FOR ALL NEURONS OSCILLATING APPRECIABLY
		if (minoutput[i] <(burstthreshold-.05)) {
			if (maxoutput[i]>burstthreshold) {
				fitness += 0.05;
				criteriamet += 1;
			}
		}
	}
		
	if (PDstartcount < 3){
		cout << "unable to find two full cycles; may want to increase transient, lengthen runtime, or speed up slowest timescale" << endl;
	}
	else{
		int PDend = 0;
		int PDendcount = 0;
		int LPstart = 0;
		int LPstartcount = 0;
		int LPend = 0;
		int PYstart = 0;
		int PYstartcount = 0;
		int PYend = 0;
		for (int step=PDstarts(1); step<=PDstarts(2); step ++){
			if (PDendcount == 0){
				if (OutputHistory[step][3]>burstthreshold){
					if (OutputHistory[step+1][3]<burstthreshold){
						PDend = step;
						PDendcount ++;
						//cout << "PDend";
					}
				}
			}
			if (LPstartcount == 0){
				if (OutputHistory[step][1]<burstthreshold){
					if (OutputHistory[step+1][1]>burstthreshold){
						LPstart = step;
						LPstartcount ++;
						//cout << "LPstart";
					}
				}
			}
			if (PYstartcount == 0){
				if (OutputHistory[step][2]<burstthreshold){
					if (OutputHistory[step+1][2]>burstthreshold){
						PYstart = step;
						PYstartcount ++;
						//cout << "PYstart";
					}
				}
			}
		}
		if (LPstartcount == 1){
			for (int step=LPstart;step<=PDstarts(3);step++){
				if (OutputHistory[step][1]>burstthreshold){
					if (OutputHistory[step+1][1]<burstthreshold){
						LPend = step;
						//cout << "LPend";
						break;
					}
				}
			}
		}
		else{cout << "LPstart not found during cycle" << endl;}

		if (PYstartcount == 1){
			for (int step=PYstart;step<=PDstarts(3);step++){
				if (OutputHistory[step][2]>burstthreshold){
					if (OutputHistory[step+1][2]<burstthreshold){
						PYend = step;
						//cout << "PYend" << endl;
						break;
					}
				}
			}
		}
		else {cout << "PYstart not found during cycle" << endl;}


		if (abs(OutputHistory[PDstarts[1]][1] - OutputHistory[PDstarts[2]][1])<tolerance){      //at the two points where PD crosses up,
			if (abs(OutputHistory[PDstarts[1]][2] - OutputHistory[PDstarts[2]][2])<tolerance){  //are the other two neurons approximately in the same place?
				// 	ORDERING CRITERIA
				if (LPstart <= PYstart){
					//cout << "order1" << endl;
					fitness += 0.05;
					criteriamet += 1;
				}
				if (LPend <= PYend){
					//cout << "order2" << endl;
					fitness += 0.05;
					criteriamet += 1;
				}
				if (PDend <= LPstart){
					//cout << "order3" << endl;
					fitness += 0.05;
					criteriamet += 1;
				}
				if (criteriamet == 6){
					//cout << LPstart << ", " << LPend << ", " << PYstart <<", " << PYend << ", " <<PDstarts[1] << ", " <<PDend <<endl;
					int period = PDstarts[2] - PDstarts[1];
					double LPfoo = LPend - LPstart; 
					double LPdutycycle = LPfoo/period; //burstduration/period
					double LPdutycyclezscore = abs(LPdutycycle - .264)/.059;
					double PYfoo = PYend-PYstart;
					double PYdutycycle = PYfoo/period; //burstduration/period
					double PYdutycyclezscore = abs(PYdutycycle - .348)/.054;
					double PDfoo = PDend-PDstarts[1];
					double PDdutycycle = PDfoo/period; //burstduration/period
					double PDdutycyclezscore = abs(PDdutycycle - .385)/.040;
					double LPbar = LPstart-PDstarts[1];
					double LPstartphase = LPbar/period; //delay/period
					double LPstartphasezscore = abs(LPstartphase - .533)/.054;
					double PYbar = PYstart-PDstarts[1];
					double PYstartphase = PYbar/period; //delay/period
					double PYstartphasezscore = abs(PYstartphase - .758)/.060;
					//cout << "Period:" << period << endl;
					//cout << LPdutycyclezscore<< ", "<<PYdutycyclezscore<<", "<<PDdutycyclezscore<<", "<<LPstartphasezscore<<", "<<PYstartphasezscore<<endl;
					double average = (LPdutycyclezscore+PYdutycyclezscore+PDdutycyclezscore+LPstartphasezscore+PYstartphasezscore)/5;
					fitness += 1/(average);
				}
			}
		}
		else{
			cout << "possible multi-periodicity" << endl;
			// cout << "LPstartcount = " << LPstartcount << " ,PYstartcount = " << PYstartcount << endl;
			// NO ORDERING POINTS FOR MULTIPERIODIC
		}
	}
	return fitness;

}

// ------------------------------------
// Recovery Fitness Function
// ------------------------------------
double HPFitnessFunction(TVector<double> &genotype, RandomState &rs){

    TVector<double> par1s(1,3);
    par1s[1] = 0.1;
    par1s[2] = 2.5;
    par1s[3] = 4.5;
    TVector<double> par2s(1,2);
    par2s[1] = -5;
    par2s[2] = -8.1;
	// cout << "par vectors: " << par1s << " " << par2s << endl;
    
    // Map genootype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);
	
	// Create the agent
    TVector<int> nullwindowsize(1,N);
    nullwindowsize.FillContents(1);
    TVector<double> nulllb(1,N);
    nulllb.FillContents(0);
    TVector<double> nullub(1,N);
    nullub.FillContents(1);
    TVector<double> nullbt(1,N);
    nullbt.FillContents(1);
    TMatrix<double> nullwt(1,N,1,N);
    nullwt.FillContents(1);
	CTRNN Agent(3, nullwindowsize, nulllb, nullub, nullbt, nullwt, 16, 16);

	// Instantiate the nervous system
	char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Agent; 

	// cout << "OG fit = " << PyloricPerformance(Agent) << endl;

    // Save the original values of theta1 and theta3
    // double theta1 = Agent.NeuronBias(1);
    // double theta3 = Agent.NeuronBias(3);

    // Instantiate the HP mechanism specified in the genome

    //Eventually, will work like this
    // int k = 1;
    // for(int i=1;i<=N;i++){
    //     Agent.SetPlasticityBoundary(i,phenotype[k]);
    // }  //....

    //For now, we will hack it (only doing theta_1, theta_3)
	// cout << phenotype << endl;
    Agent.SetNeuronBiasTimeConstant(1,phenotype[1]);
    Agent.SetNeuronBiasTimeConstant(3,phenotype[2]);
    Agent.SetPlasticityLB(1,phenotype[3]);
    Agent.SetPlasticityLB(3,phenotype[4]);
    Agent.SetPlasticityUB(1,phenotype[5]);
    Agent.SetPlasticityUB(3,phenotype[6]);
    Agent.SetSlidingWindow(1,int(phenotype[7]));
    Agent.SetSlidingWindow(3,int(phenotype[8]));

    double fitness = 0;
    for (int i=1;i<=par1s.UpperBound();i++){
        Agent.SetNeuronBias(1,par1s[i]);
        for (int j=1;j<=par2s.UpperBound();j++){
            Agent.SetNeuronBias(3,par2s[j]);
			// cout << "Biases " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << endl;
            // Initialize the outputs at 0.5 for all neurons in the circuit
            Agent.RandomizeCircuitOutput(0.5, 0.5);

            // Run the circuit for an initial transient; HP is off and fitness is not evaluated
            for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                Agent.EulerStep(StepSize,false,false);
            }

			// For testing purposes, run the circuit to see its initial pyloric fitness
			// cout << PyloricPerformance(Agent) << " ";

            // Run the circuit for a period of time with HP so the paramters can change
            for (double t = StepSize; t<= PlasticDuration; t+= StepSize){
                Agent.EulerStep(StepSize,true,true);
            }


            // Calculate the Pyloric Fitness
            fitness += PyloricPerformance(Agent);
			// cout << fitness << endl;
        }
    }  
    return fitness;
}

// ------------------------------------
// Display functions
// ------------------------------------
ofstream Evolfile;
ofstream BestIndividualsFile;

void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	
	//cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
	Evolfile << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
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
	BestIndividualsFile << s.BestPerformance() << endl << endl;

	cout << trial << "finished" << endl;

	trial ++;
}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	Evolfile.open("evol.dat");
	BestIndividualsFile.open("bestinds.dat");
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

		s.SetEvaluationFunction(HPFitnessFunction);
		s.ExecuteSearch();

		// ifstream genefile("best.gen.dat");
		// TVector<double> genotype(1, VectSize);
		// genefile >> genotype;
		
	}
	Evolfile.close();
	BestIndividualsFile.close();

  return 0;
}
