// --------------------------------------------------------------
//  Evaluate Pyloric fitness at many parameter values on a grid 
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 200; //in seconds
const double RunDuration = 200; //in seconds
const double StepSize = 0.01;
const int RunSteps = RunDuration/StepSize; // in steps

// Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .01; //for detecting double periodicity

// Parameter space resolution
const double par1min = 0.0;
const double par1max = 5.0;
const double par1step = .01;
const double par2min = -9;
const double par2max = 0;
const double par2step = .01;

// Nervous system params
const int N = 3;
// const double WR = 16.0; 
// const double BR = 16.0; //(WR*N)/2; //<-for allowing center crossing
// const double TMIN = 1; 
// const double TMAX = 2; 

// Plasticity parameters
// const int WS = 120;		// Window Size of Plastic Rule (in steps size) (so 1 is no window)
// const double B = 0.1; 		// Plasticity Low Boundary (symmetric)
// const double BT = 20.0;		// Bias Time Constant
// const double WT = 40.0;		// Weight Time Constant

int	VectSize = N*N + 2*N;

// ------------------------------------
// Pyloric-like Fitness function
// ------------------------------------
double PyloricFitnessFunction(CTRNN Agent)
{
	TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,RunSteps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);
	double fitness = 0.0;

	// Run the circuit for an initial transient; HP is off and fitness is not evaluated
	for (double t = StepSize; t <= TransientDuration; t += StepSize) {
		Agent.EulerStep(StepSize,false,false); //no HP
	}

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness while HP is turned OFF.
	int temp = 0;
	for (double time = StepSize; time <= RunDuration; time += StepSize) {
		temp += 1;
		for (int i = 1; i <= N; i += 1) {
			OutputHistory[temp][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,false,false);
	}
	int criteriamet = 0;
	for (int i = 1; i <= N; i += 1) {
		// SHORT HAND FOR ALL NEURONS OSCILLATING APPRECIABLY around the threshold of .5
		if (minoutput[i] <(burstthreshold-.05)) {
			if (maxoutput[i]>burstthreshold) {
				fitness += 0.05;
				criteriamet += 1;
			}
		}
	}
	
	if (criteriamet == 3){
		int PDstartcount = 0;
		TVector<int> PDstarts(1,3);
		PDstarts.FillContents(0);
		//LOCATE SECOND TO LAST FULL CYCLE of PD
		for (int step = RunSteps; step >= 2; step --) {
			if (PDstartcount < 3){
				if (OutputHistory[step][3] > burstthreshold){
					if (OutputHistory[step-1][3] < burstthreshold){
						PDstarts[3-PDstartcount] = step;
						PDstartcount += 1;
						//cout << "PDstarts";
					}
				}
			}
			else{
				break;
			}
		}
		if (PDstartcount < 3){
			cout << "unable to find two full cycles; may want to increase transient, lengthen runtime, or speed up slowest timescale" << endl;
		}
		else{
			int PDend = 0;
			int PDendcount = 0;
			int LPstart = 0;
			int LPstartcount = 0;
			int LPend = 0;
			int PYstart = 0;
			int PYstartcount = 0;
			int PYend = 0;
			for (int step=PDstarts(1); step<=PDstarts(2); step ++){
				if (PDendcount == 0){
					if (OutputHistory[step][3]>burstthreshold){
						if (OutputHistory[step+1][3]<burstthreshold){
							PDend = step;
							PDendcount ++;
							//cout << "PDend";
						}
					}
				}
				if (LPstartcount == 0){
					if (OutputHistory[step][1]<burstthreshold){
						if (OutputHistory[step+1][1]>burstthreshold){
							LPstart = step;
							LPstartcount ++;
							//cout << "LPstart";
						}
					}
				}
				if (PYstartcount == 0){
					if (OutputHistory[step][2]<burstthreshold){
						if (OutputHistory[step+1][2]>burstthreshold){
							PYstart = step;
							PYstartcount ++;
							//cout << "PYstart";
						}
					}
				}
			}
			if (LPstartcount == 1){
				for (int step=LPstart;step<=PDstarts(3);step++){
					if (OutputHistory[step][1]>burstthreshold){
						if (OutputHistory[step+1][1]<burstthreshold){
							LPend = step;
							//cout << "LPend";
							break;
						}
					}
				}
			}
			else{cout << "LPstart not found during cycle" << endl;}

			if (PYstartcount == 1){
				for (int step=PYstart;step<=PDstarts(3);step++){
					if (OutputHistory[step][2]>burstthreshold){
						if (OutputHistory[step+1][2]<burstthreshold){
							PYend = step;
							//cout << "PYend" << endl;
							break;
						}
					}
				}
			}
			else {cout << "PYstart not found during cycle" << endl;}


			if (abs(OutputHistory[PDstarts[1]][1] - OutputHistory[PDstarts[2]][1])<tolerance){      //at the two points where PD crosses up,
				if (abs(OutputHistory[PDstarts[1]][2] - OutputHistory[PDstarts[2]][2])<tolerance){  //are the other two neurons approximately in the same place?
					// 	ORDERING CRITERIA
					if (LPstart <= PYstart){
						//cout << "order1" << endl;
						fitness += 0.05;
						criteriamet += 1;
					}
					if (LPend <= PYend){
						//cout << "order2" << endl;
						fitness += 0.05;
						criteriamet += 1;
					}
					if (PDend <= LPstart){
						//cout << "order3" << endl;
						fitness += 0.05;
						criteriamet += 1;
					}
					if (criteriamet == 6){
						//cout << LPstart << ", " << LPend << ", " << PYstart <<", " << PYend << ", " <<PDstarts[1] << ", " <<PDend <<endl;
						int period = PDstarts[2] - PDstarts[1];
						double LPfoo = LPend - LPstart; 
						double LPdutycycle = LPfoo/period; //burstduration/period
						double LPdutycyclezscore = abs(LPdutycycle - .264)/.059;
						double PYfoo = PYend-PYstart;
						double PYdutycycle = PYfoo/period; //burstduration/period
						double PYdutycyclezscore = abs(PYdutycycle - .348)/.054;
						double PDfoo = PDend-PDstarts[1];
						double PDdutycycle = PDfoo/period; //burstduration/period
						double PDdutycyclezscore = abs(PDdutycycle - .385)/.040;
						double LPbar = LPstart-PDstarts[1];
						double LPstartphase = LPbar/period; //delay/period
						double LPstartphasezscore = abs(LPstartphase - .533)/.054;
						double PYbar = PYstart-PDstarts[1];
						double PYstartphase = PYbar/period; //delay/period
						double PYstartphasezscore = abs(PYstartphase - .758)/.060;
						//cout << "Period:" << period << endl;
						//cout << LPdutycyclezscore<< ", "<<PYdutycyclezscore<<", "<<PDdutycyclezscore<<", "<<LPstartphasezscore<<", "<<PYstartphasezscore<<endl;
						double average = (LPdutycyclezscore+PYdutycyclezscore+PDdutycyclezscore+LPstartphasezscore+PYstartphasezscore)/5;
						fitness += 1/(average);
					}
				}
			}
			else{
				cout << "possible multi-periodicity" << endl;
				// cout << "LPstartcount = " << LPstartcount << " ,PYstartcount = " << PYstartcount << endl;
				// NO ORDERING POINTS FOR MULTIPERIODIC
			}
		}
	}
	// else{
	// 	cout << "not oscillating" << endl;
	// }
	return (fitness);

}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]) 
{
	// Create file to hold data
	ofstream slicefile;
	slicefile.open("slice.dat");

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
	
	// For every pair of parameter values specified, (right now, par1=theta1, par2=theta2)
	for (double par1=par1min; par1<=par1max; par1+=par1step){
		Circuit.SetNeuronBias(1,par1);
		for (double par2=par2min; par2<=par2max; par2+=par2step){
			Circuit.SetNeuronBias(3,par2);
			Circuit.RandomizeCircuitState(0,0);
			double fitness = PyloricFitnessFunction(Circuit); // Calculate pyloric fitness
			slicefile << fitness << " "; //store in the file
		}
		slicefile << endl;
	}
		
	slicefile.close();
  return 0;
}
