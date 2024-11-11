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

int	CTRNNVectSize = N*N + 2*N;

// Pyloric Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .01; //for detecting double periodicity

// Sampling Parameters

const double TMIN = .1;
const double TMAX = 2;
const double BR = 16;
const double WR = 16;
const int num_ICs = 1000;

// Mode
const bool random_mode = false;
const bool taus_set = false;

//Filenames
char Nfname[] = "./Pete.ns";
char HPfname[] = "./converted33.dat";
char Fitnessesfname[] = "./testfitness.dat";
char ICsfname[] = "./testICs.dat";
char biastrackfname[] = "./testbiastrack.dat";
char statestrackfname[] = "./teststatestrack.dat";

const bool trackstates = false;
const int trackstatesinterval = 200; //Track neural outputs for every X trials
const bool trackparams = false;
const int trackparamsinterval = 100; //Track biases for every X trials

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

    CTRNN Circuit(3);
    // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // cout << Circuit.br;

    // Set circuit parameters (start with the given pyloric solution)
    ifstream ifs;
    ifs.open(Nfname);
    if (!ifs) {
        cerr << "File not found: " << Nfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 

    // Set of random circuit parameters to pull from
    TVector<double> genotype(1,CTRNNVectSize);
    TVector<double> phenotype(1,CTRNNVectSize);

    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    if (!HPifs) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    Circuit.SetHPPhenotype(HPifs,StepSize,true);

    cout << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityLB(2) << " " << Circuit.PlasticityLB(3) << endl;

    for (int i = 0;i<num_ICs;i++){

        long randomseed = static_cast<long>(time(NULL));
        RandomState rs(randomseed+pow(i,2));
        for (int j = 1; j <= genotype.Size(); j++)
            {genotype[j] = rs.UniformRandom(-1,1);}
        
        GenPhenMapping(genotype,phenotype);

        // FULLY RANDOM MODE (well everything but the time constants)
        if(random_mode){
            int k = 1;
            if(!taus_set){
                for(int j=1;j<=N;j++){
                    Circuit.SetNeuronTimeConstant(j,phenotype(k));
                    k++;
                }
            }
            k = N+1; //start after time constants
            for(int j=1; j<=N; j++){
                Circuit.SetNeuronBias(j,phenotype(k));
                // cout << "set a bias" << endl;
                k++;
            }

            for (int j=1; j<=N; j++){
                for (int l=1; l<=N; l++){
                    Circuit.SetConnectionWeight(j,l,phenotype(k));
                    // cout << "set a weight" << endl; 
                    k ++;   
                }
            }
        }

        // RANDOM IN HP DIMENSIONS MODE (SLICE MODE)
        else{
            int k = 1; 
            for(int j=1; j<=N; j++){
                //check for biases
                if (Circuit.plasticitypars[k]==1){
                    Circuit.SetNeuronBias(j,phenotype(k+N)); //start after time constants
                    // cout << "set a bias" << endl;
                }
                k++;
            }

            //check for weights
            for (int j=1; j<=N; j++){
                for (int l=1; l<=N; l++){
                    if (Circuit.plasticitypars[k]==1){
                        Circuit.SetConnectionWeight(j,l,phenotype(k+N)); //started after time constants
                        // cout << "set a weight" << endl;
                    }
                    k++;
                }
            }
        }


        Circuit.RandomizeCircuitState(0,0);

        // Run for transient without HP
        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0);
        }

        // Record all parameters
        ICsfile << Circuit.taus << " " << Circuit.biases << " ";
        for(int j = 1; j <= N; j ++)
        {
            for(int k=1;k<=N;k++)
            {
                ICsfile << Circuit.ConnectionWeight(j,k) << " ";
            }
        }
        ICsfile << endl;

        // Run with HP for a time
        for(double t=0;t<PlasticDuration;t+=StepSize){
            if (trackparams && (i%trackparamsinterval==0)){
                biastrack << Circuit.biases;
                for(int j = 1; j <= N; j ++){
                    for(int k=1;k<=N;k++){
                        biastrack << Circuit.ConnectionWeight(j,k) << " ";
                    }
                } 
                biastrack << endl;
            }
			if (trackstates && (i%trackstatesinterval==0)){statestrack << Circuit.NeuronOutput(1) << " " << Circuit.NeuronOutput(2) << " " << Circuit.NeuronOutput(3) << endl;}
            Circuit.EulerStep(StepSize,1);
        }
        if (trackparams && (i%trackparamsinterval==0)) {biastrack << endl;}
		if (trackstates && (i%trackstatesinterval==0)) {statestrack << endl;}

        // Record again, after HP
        ICsfile << Circuit.taus << " " << Circuit.biases << " ";
        for(int j = 1; j <= N; j ++){
            for(int k=1;k<=N;k++){
                ICsfile << Circuit.ConnectionWeight(j,k) << " ";
            }
        } 
        ICsfile << endl << endl;

        // Test for Pyloricness with HP
        double fit = PyloricPerformance(Circuit);

        fitnesses << fit << endl;

        if (trackstates){statestrack << Circuit.outputs<< " " << endl;}
    }
    fitnesses.close();
    ICsfile.close();
    biastrack.close();
	statestrack.close();
    return 0;
}