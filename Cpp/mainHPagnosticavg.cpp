// ---------------------------------------------------------
// If the thing that matters most is the HP nullclines (where
// the LC is CENTERED around the HP range), and what matters
// most to the nullclines is the average value over one cycle
// (time constants and sliding windows don't change the center,
// in theory), then just calculate the average value independently
// of HP and that will tell you, agnostically to the range,
// where HP nullclines will fall in that plane
// 
// you could even find out what the average is along the pyloric
// region boundary, and predict which combinations of 0-range HPs
// will have HPs that work

// or you could calculate the averages at the peak and try to 
// engineer the best possible HP mechanism (as i attempted to before)
// -- in what cases does this work and in what cases does it not? 
// seek to explain the unsolvable solutions
// ---------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 250; //seconds without HP, made pretty long for this analysis (considering HP is off)
const double RunDuration = 50; //seconds to look for an oscillation cycle
const double StepSize = 0.005; //small step size for accuracy
const double leaving_tolerance = 0.05; // in state space for greater accuracy
const double return_tolerance = 0.025; //less than leaving tolerance

// Parameter space resolution
const double par1min = -16.0;
const double par1max = 16.0;
const double par1step = .1;
const double par2min = -16.0;
const double par2max = 16.0;
const double par2step = .1;

// Output file
char avgsfname[] = "./HPAgnosticAverage.dat";

// Nervous system params
const int N = 3;

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create files to hold data
	ofstream avgsfile;
	avgsfile.open(avgsfname);
	
    CTRNN Circuit(N);
    char fname[] = "./pyloriccircuit.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
	// Circuit should initialize with null HP & one timestep "sliding window"
	
	// For every neuron, (even recording N2 right now, even though not considering HPs that use it)
	TVector<double> acc(1,N); //vector to store the average value of each neuron
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		Circuit.SetNeuronBias(1,par1);
		// cout << par1 << endl; //don't print if on supercomputer
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			//don't have to reset biases all the time because shouldn't actually be allowed to change...
			Circuit.SetNeuronBias(3,par2);
			acc.FillContents(0); //reset to all zeros
			Circuit.RandomizeCircuitState(0,0); //resets averaging and sliding window utilities
			// run for transient
			for (double t=StepSize;t<=TransientDuration;t+=StepSize){
				Circuit.EulerStep(StepSize,0);
			}
			// then run to keep track of the neuron values either over one full cycle or simply the run duration
			TVector<double> startstate(1,N);
			for (int i = 1; i <= N; i++){
				startstate[i] = Circuit.NeuronState(i);
			}
			// allow the neural outputs to leave and come back to the same point
			bool left = false;
			double dist = 0;
			TVector<double> avg(1,N);
			avg.FillContents(0.0);
			double t = StepSize;
			int stepnum = 0;
			while (t<= RunDuration){
				Circuit.EulerStep(StepSize,0);
				stepnum ++;
				t += StepSize;
				dist = 0;
				for (int i = 1; i <= N; i++){
					dist += pow(Circuit.NeuronState(i)-startstate[i],2);
					avg[i] += Circuit.NeuronOutput(i);
				}
				dist = sqrt(dist);

				if ((dist > leaving_tolerance)&&(!left)){
					left = true;
					// cout << "Limit cycle left" << endl;
				}

				if (left && dist < return_tolerance){
					// cout << "Limit cycle came back" << endl;
					break;
				}
			}
			for (int i = 1; i <= N; i++){
				avg[i] /= stepnum;
			}
			avgsfile << avg << endl;
		}
		avgsfile << endl;
	}	
	avgsfile.close();
  return 0;
}
