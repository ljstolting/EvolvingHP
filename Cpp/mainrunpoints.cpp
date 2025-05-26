// --------------------------------------------------------------
//  Track the parameters of many CTRNNs as they undergo Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 150; //Seconds with HP off
const double PlasticDuration = 100000; //Seconds with HP runnings
const double RunDuration = 150; //How long to test for pyloricness
// const double StepSize = 0.025;
const int RunSteps = RunDuration/StepSize;

// Nervous system params
const int N = 3;
const bool shiftedrho_tf = true;
int	CTRNNVectSize = N*N + 2*N;

// Pyloric Detection params
// const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
// const double tolerance = .01; //for detecting double periodicity

// Sampling Parameters

const double TMIN = .1;
const double TMAX = 2;
const double BRlb1 = -20;
const double BRub1 = 20;
const double BRlb3 = -20;
const double BRub3 = 20;
const double WR = 10;
const int num_ICs = 25;

// Mode
const bool random_mode = false; //randomize in other dimensions besides HP dimensions 
const bool taus_set = false; //in random mode, do we want the taus to be variable

//Filenames
char Nfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/pyloriccircuit.ns";
char HPfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/5/bestind.dat";
const bool HPfileisbestind = true;
// char HPfname[] = "./bestindtest.dat";
// char HPfname[] = "./Convenient HP Mechanisms/nullHP.dat";
// char HPfname[] = "./Specifically Evolved HP mechanisms/Pete/2D/33/bestind.dat";
// char Fitnessesfname[] = "./Convenient HP Mechanisms/Petefitbad.dat";
char Fitnessesfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/5/fit.dat";
// char ICsfname[] = "./Convenient HP Mechanisms/Peteicsbad.dat";
// char ICsfname[] = "./testics.dat";
char ICsfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/5/ics.dat";
// char biastrackfname[] = "./Convenient HP Mechanisms/Petebiastrackbad.dat";
char biastrackfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/5/biastrack.dat";
// char statestrackfname[] = "./Convenient HP Mechanisms/Petestatestrackbad.dat";
// char statestrackfname[] = "./teststatestrack.dat";
char statestrackfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/70/5/statestrack.dat";

const bool trackstates = false;
const int trackstatesinterval = 50; //Track neural outputs for every X trials
const bool trackparams = true;
const int trackparamsinterval =1; //Track biases for every X trials
const int trackingstepinterval = 5; //make the tracking files smaller by only recording every Xth step

void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Bias
	for (int i = 1; i <= 2; i++) {
		phen(k) = MapSearchParameter(gen(k), BRlb1, BRub1);
		k++;
	}
    //last neuron
    phen(k) = MapSearchParameter(gen(k), BRlb3,BRub3);
    k++;
    // cout << phen(k);
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
    Circuit.ShiftedRho(shiftedrho_tf);

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
    if(HPfileisbestind){
        Circuit.SetHPPhenotypebestind(HPifs,StepSize,true);
    }
    else{
        Circuit.SetHPPhenotype(HPifs,StepSize,true);
    }

    cout << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityLB(2) << " " << Circuit.PlasticityLB(3) << " " << Circuit.PlasticityUB(1) << " " << Circuit.PlasticityUB(2) << " " << Circuit.PlasticityUB(3) << endl;
    cout << Circuit.SlidingWindow(1) << endl;

    for (int i = 0;i<num_ICs;i++){
        long randomseed = static_cast<long>(time(NULL));
        // long randomseed = 6942069420; //if need repeats or direct compare
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
                    // cout << phenotype(k+N);
                    // Circuit.SetNeuronBias(j,-10); //if want a specific value
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


        Circuit.RandomizeCircuitOutput(0.5,0.5);
        Circuit.WindowReset();

        // Run for transient without HP
        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,false);
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
        int tstep = 0;
        for(double t=0;t<PlasticDuration;t+=StepSize){
            if (trackparams && (i%trackparamsinterval==0) && (tstep % trackingstepinterval == 0)){
                for(int j = 1; j<= Circuit.plasticitypars.Sum(); j++){
                    biastrack << Circuit.ArbDParam(j) << " "; //record only the parameters that are changing throughout the run
                }
                biastrack << endl;
            }
			if (trackstates && (i%trackstatesinterval==0) && (tstep%trackingstepinterval==0)){
                for (int j = 1; j <= Circuit.size; j++){
                    statestrack << Circuit.NeuronOutput(j) << " ";
                }
                statestrack << endl;}
            Circuit.EulerStep(StepSize,true);
            tstep ++;
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
        // double fit = PyloricPerformance(Circuit,biastrack,statestrack);
        double fit = PyloricPerformance(Circuit);

        fitnesses << fit << endl;
        // if(fit > .3){cout << "pyloric found" << endl;}
        fitnesses << Circuit.rhos << endl << endl; //proxy for whether HP is satisfied at the end, or whether it just ran into a boundary or is in a limit cycle

        if (trackstates){statestrack << Circuit.outputs<< " " << endl;}
    }
    fitnesses.close();
    ICsfile.close();
    biastrack.close();
	statestrack.close();
    return 0;
}