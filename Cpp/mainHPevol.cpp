// --------------------------------------------------------------
//  Evolve an HP mechanism to recover pyloric function best from 
//  a set of prescribed starting conditions
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 1000; //seconds without HP
const double PlasticDuration = 2000; //seconds allowing HP to act
const double TestDuration = 500; //maximum number of seconds allowed to test pyloric performance -- can be with HP still on
const bool HPtest = true;       //does HP remain on during test (shouldn't matter if platicity time constants are slow enough)
const double StepSize = 0.01;
const int TestSteps = TestDuration/StepSize; // in steps

// Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .1; //for detecting double periodicity

// EA params
const int POPSIZE = 50;
const int GENS = 100;
const int trials = 10;    // number of times to run the EA from random starting pop
const double MUTVAR = 0.1;
const double CROSSPROB = 0.0;
const double EXPECTED = 1.1;
const double ELITISM = 0.1;
const double scaling_factor = 200; // what to multiply each discrete critereon by (osc and ordering) to value them more in HP evolution

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
ofstream pyloricoutput;
// ------------------------------------
// Pyloric-like performance
// ------------------------------------
double PyloricPerformance(CTRNN &Agent)
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
	int tstep = 1;
	double t = StepSize;
	int PDstartcount = 0;
	TVector<int> PDstarts(1,3);
	PDstarts.FillContents(0);
	// cout << "PDstarts initialized"<< endl;
	while (tstep <= TestSteps && PDstartcount < 3) {
		// outputhistfile << t << " ";
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			// outputhistfile << Agent.NeuronOutput(i) << " ";
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}
		// outputhistfile << endl;

		Agent.EulerStep(StepSize,HPtest,0); // adapt only biases, not weights

		// cout << OutputHistory(tstep,3) << endl << endl;
		// cout << Agent.NeuronOutput(3) << endl;

		//Check for PD start
		if (OutputHistory(tstep,3) < burstthreshold && Agent.NeuronOutput(3) > burstthreshold){
			PDstartcount += 1;
			PDstarts[PDstartcount] = tstep;
			//cout << "PDstarts";
		}
		tstep += 1;
		t += StepSize;
	}
	// cout << "PD starts identified" << endl;
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
	if (criteriamet < 3){
		return fitness;
	}

	if (PDstartcount < 3){
		cout << "unable to find two full cycles; may want to increase transient, lengthen runtime, or speed up slowest timescale" << endl;
		return fitness;
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
				if (OutputHistory(step,3)>burstthreshold){
					if (OutputHistory(step+1,3)<burstthreshold){
						PDend = step;
						PDendcount ++;
						//cout << "PDend";
					}
				}
			}
			if (LPstartcount == 0){
				if (OutputHistory(step,1)<burstthreshold){
					if (OutputHistory(step+1,1)>burstthreshold){
						LPstart = step;
						LPstartcount ++;
						//cout << "LPstart";
					}
				}
			}
			if (PYstartcount == 0){
				if (OutputHistory(step,2)<burstthreshold){
					if (OutputHistory(step+1,2)>burstthreshold){
						PYstart = step;
						PYstartcount ++;
						//cout << "PYstart";
					}
				}
			}
		}
		if (LPstartcount == 1){
			for (int step=LPstart;step<=PDstarts(3);step++){
				if (OutputHistory(step,1)>burstthreshold){
					if (OutputHistory(step+1,1)<burstthreshold){
						LPend = step;
						//cout << "LPend";
						break;
					}
				}
			}
		}
		else{cout << "LPstart not found during cycle" << endl; return fitness;}

		if (PYstartcount == 1){
			for (int step=PYstart;step<=PDstarts(3);step++){
				if (OutputHistory(step,2)>burstthreshold){
					if (OutputHistory(step+1,2)<burstthreshold){
						PYend = step;
						//cout << "PYend" << endl;
						break;
					}
				}
			}
		}
		else {cout << "PYstart not found during cycle" << endl; return fitness;}


		if ((abs(OutputHistory(PDstarts[1],1) - OutputHistory(PDstarts[2],1))<tolerance)&&(abs(OutputHistory(PDstarts[1],2) - OutputHistory(PDstarts[2],2))<tolerance)){
			// at the two points where PD crosses up, are the other two neurons approximately in the same place?
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
				// cout << "Period:" << period << endl;
				// cout << LPstart << " " << LPend << " " << PYstart << " " << PYend << " " << PDstarts[1] << " " << PDend << endl;
				// cout << LPdutycyclezscore<< ", "<<PYdutycyclezscore<<", "<<PDdutycyclezscore<<", "<<LPstartphasezscore<<", "<<PYstartphasezscore<<endl;
				double average = (LPdutycyclezscore+PYdutycyclezscore+PDdutycyclezscore+LPstartphasezscore+PYstartphasezscore)/5;
				fitness += 1/(average);
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
	// cout << phenotype << endl;
	
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
	ifs.close();

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
	
	//////////////////////////
	// ifstream HP1;
	// HP1.open("HP1.gn");

	// ifstream HPnull;
	// HPnull.open("HPnull.gn");
	
	// Agent.SetHPPhenotype(HP1);
	// HP1.close();
	// HPnull.close();
	//////////////////////////

	// cout << "parameters" << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
    double fitness = 0;
    for (int i=1;i<=par1s.UpperBound();i++){
        Agent.SetNeuronBias(1,par1s[i]);
        for (int j=1;j<=par2s.UpperBound();j++){
            Agent.SetNeuronBias(3,par2s[j]);
			// cout << "parameters " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
            // Initialize the outputs at 0.5 for all neurons in the circuit
            Agent.RandomizeCircuitOutput(0.5, 0.5);

            // Run the circuit for an initial transient; HP is off and fitness is not evaluated
            for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                Agent.EulerStep(StepSize,false,false);
            }

			// For testing purposes, run the circuit to see its initial pyloric fitness
			// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << " " << PyloricPerformance(Agent) << " ";

            // Run the circuit for a period of time with HP so the paramters can change
            for (double t = StepSize; t<= PlasticDuration; t+= StepSize){
                Agent.EulerStep(StepSize,true,false);  //set to only adapt biases, not weights
            }


            // Calculate the Pyloric Fitness
			double fit = PyloricPerformance(Agent);
			// Transform it so Pyloricness at all is worth a lot
			if (fit <= .3){fit = fit*scaling_factor;}
			else {fit = fit+60;}
			// cout << fit << endl;
			fitness += fit;
			
            // diagnostic that determined the weights were being changed
			// cout << "parameters" << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
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
		
	}
	Evolfile.close();
	BestIndividualsFile.close();


	// TESTING CONDITION
	// RandomState rs;
	// pyloricoutput.open("pylorichist.dat");
	// double perf = HPFitnessFunction(genotype,rs);

	// Manual copy of HPfitness function
	// Load the base CTRNN parameters
    // TVector<int> Window_Sizes(1,N);
    // Window_Sizes.FillContents(1);

    // TVector<double> Lower_Bounds(1,N);
    // Lower_Bounds.FillContents(0);
    // Lower_Bounds[2] = 0;

    // TVector<double> Upper_Bounds(1,N);
    // Upper_Bounds.FillContents(1);
    // Upper_Bounds[2] = 1;

    // TVector<double> Btaus(1,N);
    // Btaus.FillContents(20);

    // TMatrix<double> Wtaus(1,N,1,N);
    // Wtaus.FillContents(40);

    // // Set HP parameters
    // CTRNN Circuit(3, Window_Sizes, Lower_Bounds, Upper_Bounds, Btaus, Wtaus, 16, 16);
    // // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // // cout << Circuit.br;
    // char fname[] = "Pete.ns";
    // ifstream ifs;
    // ifs.open(fname);
    // if (!ifs) {
    //     cerr << "File not found: " << fname << endl;
    //     exit(EXIT_FAILURE);
    // }
    // ifs >> Circuit; 

	// // ifstream HP1;
	// // HP1.open("HP1.gn");

	// // ifstream HPnull;
	// // HPnull.open("HPnull.gn");

	// // Circuit.SetHPGenome(HP1);
	// // cout << Circuit.PlasticityLB(1) << Circuit.PlasticityUB(1) << endl;

	// TVector<double> par1s(1,3);
    // par1s[1] = 0.1;
    // par1s[2] = 2.5;
    // par1s[3] = 4.5;
    // TVector<double> par2s(1,2);
    // par2s[1] = -5;
    // par2s[2] = -8.1;

	// double perf = 0.0;
    // for (int i=1;i<=par1s.UpperBound();i++){
    //     Circuit.SetNeuronBias(1,par1s[i]);
    //     for (int j=1;j<=par2s.UpperBound();j++){
    //         Circuit.SetNeuronBias(3,par2s[j]);
	// 		// cout << "Biases " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << endl;
    //         // Initialize the states at 0 for all neurons in the circuit
	// 		for (int n=1;n<=N;n++){Circuit.SetNeuronState(n,0);}

    //         // Run the circuit for an initial transient; HP is off and fitness is not evaluated
    //         for (double t = StepSize; t <= TransientDuration; t += StepSize) {
	// 			pyloricoutput << t << " ";
	// 			for(int n=1;n<=N;n++){pyloricoutput << Circuit.NeuronOutput(n)<<" ";}
	// 			pyloricoutput<<endl;
    //             Circuit.EulerStep(StepSize,false,false);
    //         }

	// 		// For testing purposes, run the circuit to see its initial pyloric fitness
	// 		// cout << PyloricPerformance(Circuit,pyloricoutput) << " ";

    //         // Run the circuit for a period of time with HP so the paramters can change
    //         for (double t = StepSize; t<= PlasticDuration; t+= StepSize){
	// 			pyloricoutput << t << " ";
	// 			for(int n=1;n<=N;n++){pyloricoutput << Circuit.NeuronOutput(n)<<" ";}
	// 			pyloricoutput<<endl;
    //             Circuit.EulerStep(StepSize,true,true);
    //         }


    //         // Calculate the Pyloric Fitness
    //         double fit = PyloricPerformance(Circuit);
	// 		// cout << fit << endl;
	// 		perf = perf + fit;
    //     }
    // }
	
	// HP1.close();
	// HPnull.close();
	
	// cout << endl << perf;
	pyloricoutput.close();
	

  return 0;
}
