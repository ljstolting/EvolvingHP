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
// const double par1min = -16;
// const double par1max = 	16;
// const double par1step = .1;
// const double par2min = -16;
// const double par2max = 16;
// const double par2step = .1;

// Nervous system params
const int N = 3;
// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create file to hold data
	ofstream slicefile;
	// slicefile.open("./Specifically Evolved HP mechanisms/Every Circuit/92/pyloricslice.dat");
	slicefile.open("./Test3DHPonPyloricSolutions/pyloricslice3D_39.dat");
	// slicefile.open("Pete_maxmindetected_HP129.dat");
	// ofstream resfile;
	// resfile.open("./Specifically Evolved HP mechanisms/Every Circuit/res_regular.dat");
	// resfile << par1min << " " << par1max << " " << par1step << endl << par2min << " " << par2max << " " << par2step << endl;
	ifstream resfile;
	resfile.open("./Test3DHPonPyloricSolutions/res.dat");
	ifstream dimsfile;
	dimsfile.open("./avgsdimensions.dat"); 

	// Load the base CTRNN parameters
    CTRNN Circuit(3);
    char fname[] = "./Specifically Evolved HP mechanisms/Every Circuit/39/pyloriccircuit.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 

	TVector<int> dims(1,N+(N*N));
	int num_dims;

	dimsfile >> dims;
	num_dims = dims.Sum();

	Circuit.ShiftedRho(true);
	Circuit.SetPlasticityPars(dims);

	TMatrix<double> resmat(1,num_dims,1,3);
	TVector<double> parvec(1,num_dims);

	for(int i=1;i<=num_dims;i++){
		for(int j=1;j<=3;j++){
			resfile >> resmat(i,j);
		}
		parvec(i) = resmat(i,1); //initialize the parameter values at the lowest given bound
	}
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta3)
	bool finished = false;
	while (!finished){
		cout << parvec << endl;
		for (int i=1;i<=num_dims;i++){
			Circuit.SetArbDParam(i,parvec(i));
		}
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
		slicefile << endl;

		
		//and then increase the value of the appropriate parameters
		parvec(num_dims)+=resmat(num_dims,3); //step the last dimension
		for (int i=(num_dims-1); i>=1; i-=1){ //start at the second to last dimension and count backwards to see if the next dimension has completed a run
			if(parvec(i+1)>resmat(i+1,2)){   //if the next dimension is over its max
				parvec(i+1) = resmat(i+1,1); //set it to its min
				parvec(i) += resmat(i,3);    //and step the current dimension
			}
		}
		if (parvec(1)>resmat(1,2)){
			finished = true;
		}
	}	
	slicefile.close();
  return 0;
}