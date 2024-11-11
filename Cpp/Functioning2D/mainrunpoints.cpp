// --------------------------------------------------------------
//  Track the parameters of many CTRNNs as they undergo Homeostatic Plasticity
// --------------------------------------------------------------
#include "../TSearch.h"
#include "CTRNN.h"
#include "../random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params
// const double TransientDuration = 1000; //Seconds with HP off
const double PlasticDuration = 10000; //Seconds with HP running

// Nervous system params
const int N = 3;

int	VectSize = N*N + 2*N;

// Pyloric Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .01; //for detecting double periodicity

// Sampling Parameters
// SLICE MODE
const double par1min = -16;
const double par1max = 16;
const double par2min = -16;
const double par2max = 16;

const double TMIN = .1;
const double TMAX = 2;
const double BR = 16;
const double WR = 16;
const int num_ICs = 1000;

//Filenames
char Nfname[] = "./TestBestonDifferentSolutions/Sven/Sven.ns";
char HPfname[] = "./HP_unevolved/HPpoint1.gn";
char Fitnessesfname[] = "./TestBestonDifferentSolutions/Sven/Fitnessespoint1.dat";
char ICsfname[] = "./TestBestonDifferentSolutions/Sven/ICspoint1.dat";
char biastrackfname[] = "./TestBestonDifferentSolutions/Sven/BiasTrackpoint1.dat";
char statestrackfname[] = "./TestBestonDifferentSolutions/Sven/StatesTrackpoint1.dat";
char detectiontrackfname[] = "./TestBestonDifferentSolutions/Sven/DetectedActivitypoint1.dat";

const bool trackstates = true;
const int trackstatesinterval = 50; //Track neural outputs for every X trials
const bool trackparams = true;
const int trackparamsinterval = 50; //Track biases for every X trials
const bool trackdetectedact = true;
const int trackdetectedactinterval = 50; //Track detected activity for every X trials

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

int main(){
    // cout << StepSize << endl;
    // Create files to hold data
	ofstream fitnesses;
    fitnesses.open(Fitnessesfname);
    ofstream ICsfile;
    ICsfile.open(ICsfname);
	ofstream biastrack;
    biastrack.open(biastrackfname);
	ofstream statestrack;
	statestrack.open(statestrackfname);
    ofstream detectedactivitytrack;
	detectedactivitytrack.open(detectiontrackfname);

    CTRNN Circuit(3);
    
    // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // cout << Circuit.br;

    // Set circuit parameters SLICE MODE
    ifstream ifs;
    ifs.open(Nfname);
    if (!ifs) {
        cerr << "File not found: " << Nfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 

    // RANDOM MODE
    // TVector<double> genotype(1,VectSize);
    // TVector<double> phenotype(1,VectSize);

    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    if (!HPifs) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    Circuit.SetHPPhenotype(HPifs,StepSize,true);

    for (int i = 1;i<=num_ICs;i++){
        // RANDOM MODE  
        // RandomState rs(pow(i,2));
        // for (int i = 1; i <= genotype.Size(); i++)
        //     {genotype[i] = rs.UniformRandom(-1,1);}
        
        // GenPhenMapping(genotype,phenotype);
        // phenotype >> Circuit;

        // SLICE MODE
        double theta1 = UniformRandom(par1min,par1max);
        double theta3 = UniformRandom(par2min,par2max);
        
        Circuit.SetNeuronBias(1,theta1);
        Circuit.SetNeuronBias(3,theta3);

        Circuit.RandomizeCircuitState(0,0);

        // Run for transient without HP
        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0,0);
        }
        // SLICE MODE
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;

        // RANDOM MODE
        // ICsfile << Circuit.taus << " " << Circuit.biases;
        // for(int j = 1; j <= N; j ++)
        // {
        //     for(int k=1;k<=N;k++)
        //     {
        //         ICsfile << Circuit.ConnectionWeight(j,k) << " ";
        //     }
        // }
        // ICsfile << endl;

        // Run with HP for a time
        for(double t=0;t<PlasticDuration;t+=StepSize){
            // cout << "t=" << t << endl;
            if (trackparams && (i%trackparamsinterval==0)){biastrack << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;}
            // if (trackparams && (i%trackparamsinterval==0)){biastrack << Circuit.biases << endl;}
			if (trackstates && (i%trackstatesinterval==0)){statestrack << Circuit.NeuronOutput(1) << " " << Circuit.NeuronOutput(2) << " " << Circuit.NeuronOutput(3) << endl;}
            if (trackdetectedact && (i%trackdetectedactinterval==0)){detectedactivitytrack << Circuit.avgoutputs << endl;}
            Circuit.EulerStep(StepSize,1,0);
        }
        if (trackparams && (i%trackparamsinterval==0)) {biastrack << endl;}
		if (trackstates && (i%trackstatesinterval==0)) {statestrack << endl;}
        if (trackdetectedact && (i%trackdetectedactinterval==0)) {detectedactivitytrack<<endl;}
        // SLICE MODE
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl << endl;

        // RANDOM MODE
        // ICsfile << Circuit.taus << " " << Circuit.biases;
        // for(int j = 1; j <= N; j ++){
        //     for(int k=1;k<=N;k++){
        //         ICsfile << Circuit.ConnectionWeight(j,k) << " ";
        //     }
        // } 
        // ICsfile << endl << endl;

        // Test for Pyloricness with HP
        double fit = PyloricPerformance(Circuit);

        fitnesses << fit << endl;

        // statestrack << Circuit.outputs<< " " << endl;
        // if ((Circuit.NeuronOutput(1) > Circuit.PlasticityLB(1)) && (Circuit.NeuronOutput(1) < Circuit.PlasticityUB(1))){
        //     statestrack << 1 << " ";
        // }
        // else {statestrack << 0 << " ";}

        // if ((Circuit.NeuronOutput(3) > Circuit.PlasticityLB(3)) && (Circuit.NeuronOutput(3) < Circuit.PlasticityUB(3))){
        //     statestrack << 1 << endl;
        // }
        // else {statestrack << 0 << endl;}
    }
    fitnesses.close();
    ICsfile.close();
    biastrack.close();
	statestrack.close();
    detectedactivitytrack.close();
    return 0;
}