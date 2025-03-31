// --------------------------------------------------------------
//  Track the parameters of a CTRNN as it undergoes Homeostatic Plasticity
//  Updated to new system 3/31/25
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

const bool debugging = true;

// Task params - can be separately defined because does not include the pyloric file
const double TransientDuration = 50; //Seconds with HP off
const double PlasticDuration = 100; //Seconds with HP running, during which states and parameters are recorded
const double StepSize = 0.01;

// Sampling Parameters
const double TMIN = .1;
const double TMAX = 2;
const double BRlb1 = -10;
const double BRub1 = 20;
const double BRlb3 = -20;
const double BRub3 = 10;
const double WR = 10;

// Nervous system params
const int N = 3;
const int CTRNNVectSize = (N*N) + (2*N);

const int randomICs = 0;    //how many randomly generated initial conditions (in the subspace of HP's action)
                            //set to 0 if only generating a timeseries for the evolved solution, without randomization
                            //all IC trajectories will be recorded

//Filenames for parameter traces
char paramoutfilename[] = "./Specifically Evolved HP mechanisms/Every Circuit/0/paramtrack.dat";
char outputsoutfilename[] = "./Specifically Evolved HP mechanisms/Every Circuit/0/outputtrack.dat";
char circuitfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/0/pyloriccircuit.ns";
char HPfname[] = "./Convenient HP Mechanisms/nullHP.dat";

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
    // Create files to hold data
	ofstream paramoutfile;
    paramoutfile.open(paramoutfilename);
    ofstream outputoutfile;
    outputoutfile.open(outputsoutfilename);
    ifstream circuitgenome;
    circuitgenome.open(circuitfname);
    if (!circuitgenome) {
        cerr << "CTRNN genome file not found: " << circuitfname << endl;
        exit(EXIT_FAILURE);
    }
    ifstream HPgenome;
    HPgenome.open(HPfname);
    if (!HPgenome) {
        cerr << "HP genome file not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }

	// Load the base CTRNN parameters and Set HP parameters
    CTRNN Circuit(3);
    circuitgenome >> Circuit; 

    // Set the proper HP parameters 
    Circuit.SetHPPhenotype(HPgenome,StepSize,true); //file-based HP setting with range encoding

    if (debugging){
        cout << "Confirm biases:" << Circuit.biases << endl;
        cout << "Confirm bounds N1: " << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityUB(1) << endl;
        cout << "Confirm bounds N3: " << Circuit.PlasticityLB(3) << " " << Circuit.PlasticityUB(3) << endl;
        cout << "Confirm SW: " << Circuit.SlidingWindow(3) << endl;
    }

    //Run the base circuit without any randomization first
    for (int i = 1; i <= N; i ++){
        Circuit.SetNeuronState(i,0);
    }

    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0);
    }

    for(double t=0;t<PlasticDuration;t+=StepSize){
        paramoutfile << Circuit.biases << endl;
        outputoutfile << Circuit.outputs << endl;
        Circuit.EulerStep(StepSize,1);
    }

    //Run designated number of randomized initial conditions, recording the action of HP
    for (int i =1; i<=randomICs; i++){
        TVector<double> genotype(1,CTRNNVectSize);
        TVector<double> phenotype(1,CTRNNVectSize);

        for (int j = 1; j <= genotype.Size(); j++)
            {genotype[j] = UniformRandom(-1,1);}
        
        GenPhenMapping(genotype,phenotype);

        int k = 1; 
        for(int j=1; j<=N; j++){
            //check for biases
            if (Circuit.plasticitypars[k]==1){
                Circuit.SetNeuronBias(j,phenotype(k+N)); //start after time constants
            }
            k++;
        }

        //check for weights
        for (int j=1; j<=N; j++){
            for (int l=1; l<=N; l++){
                if (Circuit.plasticitypars[k]==1){
                    Circuit.SetConnectionWeight(j,l,phenotype(k+N)); //started after time constants
                }
                k++;
            }
        }

        for (int j = 1; j <= N; j ++){
            Circuit.SetNeuronState(j,0);
        }

        Circuit.WindowReset();

        paramoutfile << endl;
        outputoutfile << endl;

        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0);
        }
    
        for(double t=0;t<PlasticDuration;t+=StepSize){
            paramoutfile << Circuit.biases << endl;
            outputoutfile << Circuit.outputs << endl;
            Circuit.EulerStep(StepSize,1);
        }
    }

    return 0;
}