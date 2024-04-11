// --------------------------------------------------------------
//  Track the parameters of a CTRNN as it undergoes Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 1000; //Seconds with HP off
const double PlasticDuration = 100000; //Seconds with HP running
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
char redfname[] = "./HP_unevolved/redparamtrack_handdesign.dat";
char orangefname[] = "./HP_unevolved/orangeparamtrack_handdesign.dat";
char yellowfname[] = "./HP_unevolved/yellowparamtrack_handdesign.dat";
char greenfname[] = "./HP_unevolved/greenparamtrack_handdesign.dat";
char bluefname[] = "./HP_unevolved/blueparamtrack_handdesign.dat";
char purplefname[] = "./HP_unevolved/purpleparamtrack_handdesign.dat";

char HPfname[] = "./HP_unevolved/HPhanddesign.gn";

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

	// Load the base CTRNN parameters and Set HP parameters
    CTRNN Circuit(3);
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
    Circuit.SetHPPhenotype(HPifs,StepSize);

    cout << "Confirm B1: " << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityUB(1) << endl;
    cout << "Confirm B3: " << Circuit.PlasticityLB(3) << " " << Circuit.PlasticityUB(3) << endl;
    cout << "Confirm SW: " << Circuit.SlidingWindow(3) << endl;


    // Run the red point
    Circuit.SetNeuronBias(1,.1);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
        // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    }
    // Circuit.PrintMaxMinAvgs();
    for(double t=0;t<PlasticDuration;t+=StepSize){
        redfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
        // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end red" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;

    // Run the orange point
    Circuit.SetNeuronBias(1,2.5);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<PlasticDuration;t+=StepSize){
        orangefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end orange" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;
    
    // Run the yellow point
    Circuit.SetNeuronBias(1,4.5);
    Circuit.SetNeuronBias(3,-5);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<PlasticDuration;t+=StepSize){
        yellowfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end yellow" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;
    
    // Run the green point
    Circuit.SetNeuronBias(1,4.5);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<PlasticDuration;t+=StepSize){
        greenfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end green" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;
    
    // Run the blue point
    Circuit.SetNeuronBias(1,2.5);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<PlasticDuration;t+=StepSize){
        bluefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end blue" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;

    // Run the purple point
    Circuit.SetNeuronBias(1,.1);
    Circuit.SetNeuronBias(3,-8.1);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        Circuit.EulerStep(StepSize,0,0);
    }
    for(double t=0;t<PlasticDuration;t+=StepSize){
        purplefile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    Circuit.PrintMaxMinAvgs();
    cout << "confirm end purple" << Circuit.biases(1) << " " << Circuit.biases(3) << endl;
    
    return 0;
}