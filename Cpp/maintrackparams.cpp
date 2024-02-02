// --------------------------------------------------------------
//  Track the parameters of a CTRNN as it undergoes Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 250; //Seconds with HP off
const double RunDuration = 1000; //Seconds with HP running
const double StepSize = 0.01;

// Nervous system params
const int N = 3;
const double TMIN = .1; 
const double TMAX = 2; 

// Plasticity parameters
const int WS = 0;		// Window Size of Plastic Rule (in steps size) (so 1 is no window)
const double B = 0.25; 		// Plasticity Low Boundary (symmetric)
const double BT = 20.0;		// Bias Time Constant
const double WT = 40.0;		// Weight Time Constant
const double WR = 100;      // Range that weights cannot exceed (make large to not matter)
const double BR = 100;      // Range that biases cannot exceed (make large to not matter)

int	VectSize = N*N + 2*N;

//Filenames for parameter traces
char redfname[] = "redparamtrack2.dat";
char orangefname[] = "orangeparamtrack2.dat";
char yellowfname[] = "yellowparamtrack2.dat";
char greenfname[] = "greenparamtrack2.dat";
char bluefname[] = "blueparamtrack2.dat";
char purplefname[] = "purpleparamtrack2.dat";

char HPfname[] = "HP2.gn";

// Will eventually need a function to map genotype of HP into phenotype of values

int main(){
    // Create files to hold data
	ofstream redfile;
	redfile.open(redfname);
    ofstream orangefile;
	orangefile.open(orangefname);
    ofstream yellowfile;
	yellowfile.open(yellowfname);
    ofstream greenfile;
	greenfile.open(greenfname);
    ofstream bluefile;
	bluefile.open(bluefname);
    ofstream purplefile;
	purplefile.open(purplefname);

	// Load the base CTRNN parameters
    TVector<int> Window_Sizes(1,N);
    Window_Sizes.FillContents(1);

    TVector<double> Lower_Bounds(1,N);
    Lower_Bounds.FillContents(B);
    Lower_Bounds[2] = 0;

    TVector<double> Upper_Bounds(1,N);
    Upper_Bounds.FillContents(1-B);
    Upper_Bounds[2] = 1;

    TVector<double> Btaus(1,N);
    Btaus.FillContents(BT);

    TMatrix<double> Wtaus(1,N,1,N);
    Wtaus.FillContents(WT);

    // Set HP parameters
    CTRNN Circuit(3, Window_Sizes, Lower_Bounds, Upper_Bounds, Btaus, Wtaus, BR, WR);
    // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // cout << Circuit.br;
    char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    Circuit.SetHPGenome(HPifs);

    cout << "Confirm LB1: " << Circuit.PlasticityLB(1) << endl;


    // Run the red point
    Circuit.SetNeuronBias(1,.1);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
        // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        redfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
        // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    }

    // Run the orange point
    Circuit.SetNeuronBias(1,2.5);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        orangefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    
    // Run the yellow point
    Circuit.SetNeuronBias(1,4.5);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        yellowfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    
    // Run the green point
    Circuit.SetNeuronBias(1,4.5);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        greenfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    
    // Run the blue point
    Circuit.SetNeuronBias(1,2.5);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        bluefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }

    // Run the purple point
    Circuit.SetNeuronBias(1,.1);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<RunDuration;t+=StepSize){
        purplefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    
    return 0;
}