// --------------------------------------------------------------
//  Track the parameters of many CTRNNs as they undergo Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params
// const double TransientDuration = 1000; //Seconds with HP off
const double PlasticDuration = 10000; //Seconds with HP running
const double RunDuration = 150; //How long to test for pyloricness
// const double StepSize = 0.025;
const int RunSteps = RunDuration/StepSize;

// Nervous system params
const int N = 3;

int	VectSize = N*N + 2*N;

// Pyloric Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .01; //for detecting double periodicity

// Sampling Parameters
const double par1min = -16;
const double par1max = 16;
const double par2min = -16;
const double par2max = 16;
const int num_ICs = 1000;

//Filenames
char Nfname[] = "Ronn.ns";
char HPfname[] = "./33/bestind.dat";
char Fitnessesfname[] = "RonnFinalFitnesses33.dat";
char ICsfname[] = "RonnICs33.dat";
char biastrackfname[] = "RonnBiasTrack33.dat";
char statestrackfname[] = "RonnStatesTrack33.dat";

const bool trackstates = false;
const int trackstatesinterval = 200; //Track neural outputs for every X trials
const bool trackparams = true;
const int trackparamsinterval = 20; //Track biases for every X trials

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

    CTRNN Circuit(3);
    // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // cout << Circuit.br;
    ifstream ifs;
    ifs.open(Nfname);
    if (!ifs) {
        cerr << "File not found: " << Nfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    if (!HPifs) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    Circuit.SetHPPhenotypebestind(HPifs,StepSize,true);

    for (int i = 1;i<=num_ICs;i++){
        double theta1 = UniformRandom(par1min,par1max);
        double theta3 = UniformRandom(par2min,par2max);
        
        Circuit.SetNeuronBias(1,theta1);
        Circuit.SetNeuronBias(3,theta3);
        Circuit.RandomizeCircuitState(0,0);

        // Run for transient without HP
        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0,0);
        }
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        // Run with HP for a time
        for(double t=0;t<PlasticDuration;t+=StepSize){
            if (trackparams && (i%trackparamsinterval==0)){biastrack << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;}
			if (trackstates && (i%trackstatesinterval==0)){statestrack << Circuit.NeuronOutput(1) << " " << Circuit.NeuronOutput(2) << " " << Circuit.NeuronOutput(3) << endl;}
            Circuit.EulerStep(StepSize,1,0);
        }
        if (trackparams && (i%trackparamsinterval==0)) {biastrack << endl;}
		if (trackstates && (i%trackstatesinterval==0)) {statestrack << endl;}
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl << endl;

        // Test for Pyloricness with HP
        double fit = PyloricPerformance(Circuit);

        fitnesses << fit << endl;
    }
    fitnesses.close();
    ICsfile.close();
    biastrack.close();
	statestrack.close();
    return 0;
}