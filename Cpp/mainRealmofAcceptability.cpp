// --------------------------------------------------------------
//  Evaluate the "Realm of Acceptability" (where parameters are not changing) for a given HP mechanism
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 50; //seconds without HP
const double RunDuration = 50; //seconds to test whether HP activates
const double StepSize = 0.01;

// Parameter space resolution
const double par1min = 0.0;
const double par1max = 10.0;
const double par1step = .05;
const double par2min = -15.;
const double par2max = 0.0;
const double par2step = .05;

// Run mode
const bool accelerated = false;  //stops evaluation if HP has activated for both neurons -- MESSES UP AVERAGE DETECTED

// HP genome file
char HPname[] = "./HP_unevolved/HPtest.gn"; 

// Nervous system params
const int N = 3;


// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create files to hold data
	ofstream realmacceptabilityfile;
	realmacceptabilityfile.open("realmacceptability.dat");
	ofstream maxmindetectedfile;
	maxmindetectedfile.open("maxmindetected.dat");

	// Load the base CTRNN parameters (no HP)
	TVector<int> nullwindowsize(1,N);
    nullwindowsize.FillContents(1);
    TVector<double> nulllb(1,N);
    nulllb.FillContents(0);
    TVector<double> nullub(1,N);
    nullub.FillContents(1);
    TVector<double> nullbt(1,N);
    nullbt.FillContents(1);
    TMatrix<double> nullwt(1,N,1,N);
    nullwt.FillContents(1);
	
    CTRNN Circuit(3, nullwindowsize, nulllb, nullub, nullbt, nullwt, 16, 16);
    char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
	ifstream HPphen;
	HPphen.open(HPname);
	if (!HPphen) {
        cerr << "File not found: " << HPname << endl;
        exit(EXIT_FAILURE);
    }

	Circuit.SetHPPhenotype(HPphen,StepSize);
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta3)
	TVector<bool> acc(1,2); //vector to store a boolean of whether HP was activated for each parameter
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			acc.FillContents(0); //reset to all false
			Circuit.SetMaxavg(1,0); //reset max and min detected
			Circuit.SetMinavg(1,1);
			Circuit.SetMaxavg(3,0);
			Circuit.SetMinavg(3,1);
			Circuit.SetNeuronBias(1,par1);
			Circuit.SetNeuronBias(3,par2);
			Circuit.RandomizeCircuitState(0,0);
			// run for a short transient
			for (double t=StepSize;t<=TransientDuration;t+=StepSize){
				Circuit.EulerStep(StepSize,0,0);
			}
			// then run to see whether HP would ever turn on
			// -- NOT the same as running *with* HP, so it is still directly analagous to the static map
			// -- OR IS IT? if HP does turn on, at any point after the tranisent, then it both becomes a different system and also is already flagged. OH well
			for (double t=StepSize;t<=RunDuration;t+=StepSize){
				if (acc(1)==1 && acc(2)==1 && accelerated){break;}
				Circuit.EulerStep(StepSize,0,0);
				Circuit.RhoCalc();
				//Specific to theta1/theta3 changing
				if (Circuit.NeuronRho(1)!=0){
					acc(1)=1;
				} 
				if (Circuit.NeuronRho(3)!=0){
					acc(2)=1;
				}
			}
			realmacceptabilityfile << acc << endl; //store in the file
			maxmindetectedfile << Circuit.maxavg << endl;
			maxmindetectedfile << Circuit.minavg << endl << endl;
		}
		realmacceptabilityfile << endl;
	}	
	realmacceptabilityfile.close();
  return 0;
}
