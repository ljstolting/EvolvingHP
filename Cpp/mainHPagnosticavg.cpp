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

// UPDATE 1/23/25 (post ICMNS submission)
// Averaging is only informative for HP behavior @ fixed points.
// For limit cycles, I've devised a new HP-agnostic measure which is based on the 
// rhythm's duty cycle, in combination with its average activity
// level, along with the assumption that oscillations mostly pause
// in the very top and very bottom of the output range (due to sigmoid
// nonlinearity). The hope is that this HP-agnostic measure will
// be able to predict the 
// ---------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//#define PRINTOFILE

// Task params
// const double TransientDuration = 250; //seconds without HP, made pretty long for this analysis (considering HP is off)
const double RunDuration = 50; //seconds to look for an oscillation cycle
// const double StepSize = 0.005; //small step size for accuracy
const double leaving_tolerance = 0.05; // in state space for greater accuracy
const double return_tolerance = 0.025; //less than leaving tolerance

// Parameter space resolution
const double par1min = 0.0;
const double par1max = 16.0;
const double par1step = .05;
const double par2min = 0.0;
const double par2max = 16.0;
const double par2step = .05;

// HP mode
const bool rhoshifted = true; //shifted rho is the new version of HP, where slopes are not adjusted
							   //the shifted rho HP is the only version for which I've come up with a
							   // potentially valid HP-agnostic approximation

// Output file
char avgsfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/18/HPAgnosticAverage_highres_newrho.dat";
char avgs_multistabilityfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/18/HPAgnosticAverage_highres_multistabilitytest.dat";

// Nervous system params
const int N = 3;

// Multistability detection mode: to explain trajectory divergences that occur due to multiperiodicity
// without actually finding the limit set at every point, we simply run the average detection mechanism
// from multiple starting neural states. In the base case, we will run them from a square grid 
// FUTURE MULTISTABILITY MODE SHOULD START FROM THE ENDPOINT OF THE LAST CIRCUIT (i.e. just don't reset it)
const bool multistable_mode = false;


// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create files to hold data
	ofstream avgsfile;
	avgsfile.open(avgsfname);
	ofstream avgsfilemultistability;
	avgsfilemultistability.open(avgs_multistabilityfname);
	
    CTRNN Circuit(N);
    char fname[] = "./Specifically Evolved HP mechanisms/Every Circuit/18/pyloriccircuit.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
	Circuit.ShiftedRho(rhoshifted);

	int num_pts;
	TMatrix<double> ptlist;
	if (multistable_mode){
		int resolution = 2;
		TVector<double> output_vals(1,resolution);
		output_vals[1] = .25;
		output_vals[2] = .5;
		num_pts = pow(resolution,N);
		ptlist.SetBounds(1,num_pts,1,N);
		PointGrid(ptlist,output_vals);
	}

	else{
		num_pts = 1;
		ptlist.SetBounds(1,num_pts,1,N);
		ptlist.FillContents(.5);
		// ptlist.FillContents(0); //just for now, setting circuit state to zero to align with pyloric plane
	}

	TMatrix<double> unique_endpoints(1,2,1,N);
	
	// For every neuron, (even recording N2 right now, even though not considering HPs that use it)
	TVector<double> acc(1,N); //vector to store the value of the proxy expression for each neuron
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		cout << par1 << endl;
		Circuit.SetNeuronBias(1,par1);
		// cout << par1 << endl; //don't print if on supercomputer
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			//don't have to reset biases all the time because shouldn't be allowed to change...
			Circuit.SetNeuronBias(3,par2);
			unique_endpoints.FillContents(0);
			int num_unique_endpoints = 1;
			acc.FillContents(0); //reset to all zeros
			for (int ic=1;ic<=num_pts;ic++){
				for (int neuron=1;neuron<=N;neuron++){
					Circuit.SetNeuronOutput(neuron,ptlist(ic,neuron)); //just for now, setting circuit state to align with pyloric fitness plane
					// Circuit.SetNeuronState(neuron,ptlist(ic,neuron)); 
					Circuit.WindowReset();
				}
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
				TVector<double> dutycycle(1,N);
				dutycycle.FillContents(0.0);
				TVector<bool> belowthresh(1,N);
				belowthresh.FillContents(0);
				while (t<= RunDuration){
					Circuit.EulerStep(StepSize,0);
					stepnum ++;
					t += StepSize;
					dist = 0;
					for (int i = 1; i <= N; i++){
						dist += pow(Circuit.NeuronState(i)-startstate[i],2);
						avg[i] += Circuit.NeuronOutput(i);
						if (Circuit.NeuronOutput(i) > 0.5){
							dutycycle[i] += 1;
						}
						else{
							if (!belowthresh(i)){
								belowthresh(i) = true;
							}
						}
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
					// make into true average and duty cycle by dividing by the number of steps in the cycle
					avg[i] /= stepnum;
					dutycycle[i] /= stepnum;
				}

				if (ic == num_pts){
					avgsfile << avg << endl;
					for (int i = 1; i <= N; i++){
						if ((dutycycle(i) > 0) && belowthresh(i)){ //if the neuron went above and below the threshold
							avgsfile << dutycycle(i) << " ";
						}
						else{
							avgsfile << 0 << " ";
						}
					}
					avgsfile << endl << endl;

				}

				bool uniqueness_flag = true;
				for (int i = 1; i < num_unique_endpoints; i++){
					// cout << avg << endl << unique_endpoints(i,1) << " " << unique_endpoints(i,3) << endl;
					bool points_equal = true;
					for (int j = 1; j <= N; j++){
						if (abs(avg[j] - unique_endpoints(i,j))>return_tolerance){
							points_equal = false;
						}
					}
					// cout << points_equal << endl << endl;
					if (points_equal == true){uniqueness_flag = false;}
				}

				if (uniqueness_flag && num_unique_endpoints <= 2){ //pretty sure you can only have two stable solutions in a 3N CTRNN...
					for (int i = 1; i <= N; i++){
						unique_endpoints(num_unique_endpoints,i) = avg[i];
					}
					num_unique_endpoints ++;
				}

			}
			for (int i = unique_endpoints.RowLowerBound(); i <= unique_endpoints.RowUpperBound(); i++){
				for (int j = 1; j <= N; j++){
					if (multistable_mode) {avgsfilemultistability << unique_endpoints(i,j) << " ";}
				}
				if (multistable_mode) {avgsfilemultistability << endl;}
			}
			if (multistable_mode){avgsfilemultistability << endl;}
		}
		if(multistable_mode){avgsfilemultistability << endl;}
	}	
	avgsfile.close();
	avgsfilemultistability.close();
  return 0;
}
