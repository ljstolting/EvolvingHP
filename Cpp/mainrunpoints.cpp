// --------------------------------------------------------------
//  Track the parameters of many CTRNNs as they undergo Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 1000; //Seconds with HP off
const double PlasticDuration = 10000; //Seconds with HP running
const double RunDuration = 100; //How long to test for pyloricness
const double StepSize = 0.025;
const int RunSteps = RunDuration/StepSize;

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

// Pyloric Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .01; //for detecting double periodicity

// Sampling Parameters
const double par1min = -16;
const double par1max = 16;
const double par2min = -16;
const double par2max = 16;
const int num_ICs = 1000;

//Filenames
char Nfname[] = "Sven.ns";
char HPfname[] = "./HP_unevolved/HPhanddesignstat.gn";
char Fitnessesfname[] = "SvenFinalFitnesseslong.dat";
char ICsfname[] = "SvenICslong.dat";
char biastrackfname[] = "SvenBiasTracklong.dat";
char statestrackfname[] = "SvenStatesTracklong.dat";

const bool trackstates = false;
const int trackstatesinterval = 200; //Track neural outputs for every X trials
const bool trackparams = false;
const int trackparamsinterval = 20; //Track biases for every X trials

// ------------------------------------
// Pyloric-like Fitness function -- 
// Version where you just run it from where you are (already equilibrated) with HP on
// ------------------------------------
double PyloricFitnessFunction(CTRNN Agent, bool HPtest)
{
	TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,RunSteps,1,N);
	OutputHistory.FillContents(0.0);
	double fitness = 0.0;

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness while HP is on or off depending on input
	int temp = 0;
	for (double time = StepSize; time <= RunDuration; time += StepSize) {
		temp += 1;
		for (int i = 1; i <= N; i += 1) {
			OutputHistory[temp][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,HPtest,false);
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
	return (fitness);

}

int main(){
    // Create files to hold data
	ofstream fitnesses;
    fitnesses.open(Fitnessesfname);
    ofstream ICsfile;
    ICsfile.open(ICsfname);
	ofstream biastrack;
    biastrack.open(biastrackfname);
	ofstream statestrack;
	statestrack.open(statestrackfname);

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
    ifstream ifs;
    ifs.open(Nfname);
    if (!ifs) {
        cerr << "File not found: " << Nfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    if (!HPifs) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    Circuit.SetHPPhenotype(HPifs,StepSize);

    for (int i = 1;i<=num_ICs;i++){
        double theta1 = UniformRandom(par1min,par1max);
        double theta3 = UniformRandom(par2min,par2max);
        
        Circuit.SetNeuronBias(1,theta1);
        Circuit.SetNeuronBias(3,theta3);
        Circuit.RandomizeCircuitState(0,0);

        // Run for transient without HP
        for(double t=0;t<TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0,0);
        }
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        // Run with HP for a time
        for(double t=0;t<PlasticDuration;t+=StepSize){
            if (trackparams && (i%trackparamsinterval==0)){biastrack << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;}
			if (trackstates && (i%trackstatesinterval==0)){statestrack << Circuit.NeuronOutput(1) << " " << Circuit.NeuronOutput(2) << " " << Circuit.NeuronOutput(3) << endl;}
            Circuit.EulerStep(StepSize,1,0);
        }
        if (trackparams && (i%trackparamsinterval==0)) {biastrack << endl;}
		if (trackstates && (i%trackstatesinterval==0)) {statestrack << endl;}
        ICsfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl << endl;

        // Test for Pyloricness with HP
        double fit = PyloricFitnessFunction(Circuit,true);

        fitnesses << fit << endl;
    }
    fitnesses.close();
    ICsfile.close();
    biastrack.close();
	statestrack.close();
    return 0;
}