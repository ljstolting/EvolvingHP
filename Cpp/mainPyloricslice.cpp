// --------------------------------------------------------------
//  Evaluate Pyloric fitness at many parameter values on a grid 
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params -- defined in pyloric.h file
const double TransientDuration = 100; //in seconds
// const double RunDuration = 250; //in seconds
// const int RunSteps = RunDuration/StepSize; // in steps

// Parameter space resolution
const double par1min = -16;
const double par1max = 	16;
const double par1step = .1;
const double par2min = -16;
const double par2max = 16;
const double par2step = .1;

// Nervous system params
const int N = 3;
// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create file to hold data
	ofstream slicefile;
	slicefile.open("./Specifically Evolved HP mechanisms/Every Circuit/92/pyloricslice.dat");
	// slicefile.open("Pete_maxmindetected_HP129.dat");
	ofstream resfile;
	resfile.open("./Specifically Evolved HP mechanisms/Every Circuit/res_regular.dat");
	resfile << par1min << " " << par1max << " " << par1step << endl << par2min << " " << par2max << " " << par2step << endl;

	// Load the base CTRNN parameters
    CTRNN Circuit(3);
    char fname[] = "./Specifically Evolved HP mechanisms/Every Circuit/92/pyloriccircuit.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
	Circuit.ShiftedRho(true);
	// cout << "2" << Circuit.biases << endl;
	// cout << Circuit.l_boundary << " " << Circuit.u_boundary << " " << Circuit.windowsize << endl;
	// HP should be null
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta3)
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		// cout << par1 << endl; //don't print if on supercomputer
		Circuit.SetNeuronBias(1,par1);
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			Circuit.SetNeuronBias(3,par2);
			// Circuit.RandomizeCircuitState(0,0); // resets sliding window calculation utilities, as well
			for (int neuron=1; neuron <= N; neuron ++){
				Circuit.SetNeuronOutput(neuron,.5); //puts in line with the way the average/proxy is calculated most cleanly
			}
			// cout << "3" << Circuit.biases << endl;
			// cout << Circuit.windowsize << endl << Circuit.minavg << endl << Circuit.maxavg << endl << endl;

			for (double t = StepSize; t<= TransientDuration; t+=StepSize){
				Circuit.EulerStep(StepSize,false);
			}
			// cout << "4" << Circuit.biases << endl;
			// for (double t = StepSize; t<= TestDuration; t+=StepSize){
			// 	Circuit.EulerStepAvgsnoHP(StepSize);
			// }
			// slicefile << Circuit.minavg << endl << Circuit.maxavg << endl << endl;
			slicefile << PyloricPerformance(Circuit) << " ";
			// cout << "5" << Circuit.biases << endl;
		}
		slicefile << endl;
	}	
	slicefile.close();
  return 0;
}