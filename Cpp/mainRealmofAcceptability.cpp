// ---------------------------------------------------------
// Evaluate the net parameter change that would be incurred 
// over 500 seconds (= 50 * max SW) for 
// ---------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 250; //seconds without HP
const double RunDuration = 500; //seconds to accumulate movement
const double StepSize = 0.01;

// Parameter space resolution
const double par1min = -16.0;
const double par1max = 16.0;
const double par1step = .1;
const double par2min = -16.0;
const double par2max = 16.0;
const double par2step = .1;

// HP genome file
char HPfname[] = "./33/bestind.dat"; 
char netchangefname[] = "Ronnnetchangelong.dat";

// Nervous system params
const int N = 3;

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create files to hold data
	ofstream netchangefile;
	netchangefile.open(netchangefname);
	
    CTRNN Circuit(N);
    char fname[] = "Ronn.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
	ifstream HPphen;
	HPphen.open(HPfname);
	if (!HPphen) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
	bool range_encoding = false;
	Circuit.SetHPPhenotypebestind(HPphen,StepSize,range_encoding);
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta3)
	TVector<double> acc(1,2); //vector to store the magnitude of change accumulated for each parameter during 
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		Circuit.SetNeuronBias(1,par1);
		cout << par1 << endl;
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			//don't have to reset neuron1 all the time because shouldn't actually be allowed to change...
			Circuit.SetNeuronBias(3,par2);
			acc.FillContents(0); //reset to all zeros
			Circuit.RandomizeCircuitState(0,0); //resets averaging and sliding window utilities
			// run for a short transient
			for (double t=StepSize;t<=TransientDuration;t+=StepSize){
				Circuit.EulerStep(StepSize,0,0);
			}
			// then run to keep track of the bias changes, but not actually apply them
			for (double t=StepSize;t<=RunDuration;t+=StepSize){
				Circuit.EulerStepAvgsnoHP(StepSize);
				// EVENTUALLY WHEN ALL THREE BIASES
				// for (int i = 1; i <= acc.UpperBound();i++){
				// 	acc(i) += StepSize * Circuit.RtausBiases[i] * Circuit.rhos[i];
				// }
				acc(1) += StepSize * Circuit.RtausBiases[1] * Circuit.rhos[1];
				acc(2) += StepSize * Circuit.RtausBiases[3] * Circuit.rhos[3];
			}
			netchangefile << acc << endl;
		}
		netchangefile << endl;
	}	
	netchangefile.close();
  return 0;
}
