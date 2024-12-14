// --------------------------------------------------------------
//  Given a slice data file (HP realm of acceptability), identify 
//  the points that make up the nullclines for N1 and N3, whether
//  each point is stable or unstable depending on the direction of
//  the zero-crossing, and where the nullclines cross, if anywhere. 
//  
//  Assume a 0-range HP, whereby there will not be a large swath
//  of points where no change occurs. 
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params -- defined in pyloric.h file
// const double TransientDuration = 500; //in seconds
// const double RunDuration = 250; //in seconds
// const int RunSteps = RunDuration/StepSize; // in steps

// Parameter space resolution
const double par1min = -16;
const double par1max = 16;
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
	slicefile.open("Trickster_slice.dat");
	// slicefile.open("Pete_maxmindetected_HP129.dat");

	// Load the base CTRNN parameters
    CTRNN Circuit(3);
    char fname[] = "../Pyloric CTRNN Genomes/Trickster.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 

	// Load in the HP mechanism for the purpose of specifying the sliding windows
	ifstream HPfile;
	char HPfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/0/0/bestind.dat";
	HPfile.open(HPfname);
	bool range_encoding = true;

	Circuit.SetHPPhenotype(HPfile,StepSize,range_encoding);
	// cout << Circuit.windowsize << endl;
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta3)
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		cout << par1 << endl;
		Circuit.SetNeuronBias(1,par1);
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			Circuit.SetNeuronBias(3,par2);
			Circuit.RandomizeCircuitState(0,0); // resets sliding window calculation utilities, as well
			// cout << Circuit.windowsize << endl << Circuit.minavg << endl << Circuit.maxavg << endl << endl;

			for (double t = StepSize; t<= TransientDuration; t+=StepSize){
				Circuit.EulerStep(StepSize,false);
			}
			// for (double t = StepSize; t<= TestDuration; t+=StepSize){
			// 	Circuit.EulerStepAvgsnoHP(StepSize);
			// }
			// slicefile << Circuit.minavg << endl << Circuit.maxavg << endl << endl;
			slicefile << PyloricPerformance(Circuit) << " ";
		}
		slicefile << endl;
	}	
	slicefile.close();
  return 0;
}
