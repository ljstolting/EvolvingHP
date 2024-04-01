#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

const double TransientDuration = 5000; //in seconds
const double RunDuration = 50000; //in seconds
const double StepSize = 0.01;
const int RunSteps = RunDuration/StepSize; // in steps
const double TestDuration = 500;
const int TestSteps = TestDuration/StepSize; // in steps

// Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .1; //for detecting double periodicity

const int N = 3;
int	VectSize = N*N + 2*N;

double PyloricPerformance(CTRNN &Agent)
{
    TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);
	double fitness = 0.0;

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness -- HP is either left on or turned OFF depending on specification.

	// Run the circuit until you identify 3 PD starts (capping 2 full cycles), keeping track of whether each other neuron crossed the threshold or not
	int tstep = 1;
	double t = StepSize;
	int PDstartcount = 0;
	TVector<int> PDstarts(1,3);
	PDstarts.FillContents(0);
	// cout << "PDstarts initialized"<< endl;
	while (tstep <= TestSteps && PDstartcount < 3) {
		// outputhistfile << t << " ";
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			// outputhistfile << Agent.NeuronOutput(i) << " ";
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}
		// outputhistfile << endl;

		Agent.EulerStep(StepSize,0,0); // adapt neither biases, nor weights

		// cout << OutputHistory(tstep,3) << endl << endl;
		// cout << Agent.NeuronOutput(3) << endl;

		//Check for PD start
		if (OutputHistory(tstep,3) < burstthreshold && Agent.NeuronOutput(3) > burstthreshold){
			PDstartcount += 1;
			PDstarts[PDstartcount] = tstep;
			//cout << "PDstarts";
		}
		tstep += 1;
		t += StepSize;
	}
	// cout << "PD starts identified" << endl;
	int criteriamet = 0;
	for (int i = 1; i <= N; i += 1) {
		// SHORT HAND FOR ALL NEURONS OSCILLATING APPRECIABLY
		if (minoutput[i] <(burstthreshold-.05)) {
			if (maxoutput[i]>burstthreshold) {
				fitness += 0.05;
				criteriamet += 1;
			}
		}
	}
	if (criteriamet < 3){
		return fitness;
	}

	if (PDstartcount < 3){
		cout << "unable to find two full cycles; may want to increase transient, lengthen runtime, or speed up slowest timescale" << endl;
		return fitness;
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
				if (OutputHistory(step,3)>burstthreshold){
					if (OutputHistory(step+1,3)<burstthreshold){
						PDend = step;
						PDendcount ++;
						//cout << "PDend";
					}
				}
			}
			if (LPstartcount == 0){
				if (OutputHistory(step,1)<burstthreshold){
					if (OutputHistory(step+1,1)>burstthreshold){
						LPstart = step;
						LPstartcount ++;
						//cout << "LPstart";
					}
				}
			}
			if (PYstartcount == 0){
				if (OutputHistory(step,2)<burstthreshold){
					if (OutputHistory(step+1,2)>burstthreshold){
						PYstart = step;
						PYstartcount ++;
						//cout << "PYstart";
					}
				}
			}
		}
		if (LPstartcount == 1){
			for (int step=LPstart;step<=PDstarts(3);step++){
				if (OutputHistory(step,1)>burstthreshold){
					if (OutputHistory(step+1,1)<burstthreshold){
						LPend = step;
						//cout << "LPend";
						break;
					}
				}
			}
		}
		else{cout << "LPstart not found during cycle" << endl; return fitness;}

		if (PYstartcount == 1){
			for (int step=PYstart;step<=PDstarts(3);step++){
				if (OutputHistory(step,2)>burstthreshold){
					if (OutputHistory(step+1,2)<burstthreshold){
						PYend = step;
						//cout << "PYend" << endl;
						break;
					}
				}
			}
		}
		else {cout << "PYstart not found during cycle" << endl; return fitness;}


		if ((abs(OutputHistory(PDstarts[1],1) - OutputHistory(PDstarts[2],1))<tolerance)&&(abs(OutputHistory(PDstarts[1],2) - OutputHistory(PDstarts[2],2))<tolerance)){
			// at the two points where PD crosses up, are the other two neurons approximately in the same place?
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
				// cout << "Period:" << period << endl;
				// cout << LPstart << " " << LPend << " " << PYstart << " " << PYend << " " << PDstarts[1] << " " << PDend << endl;
				// cout << LPdutycyclezscore<< ", "<<PYdutycyclezscore<<", "<<PDdutycyclezscore<<", "<<LPstartphasezscore<<", "<<PYstartphasezscore<<endl;
				double average = (LPdutycyclezscore+PYdutycyclezscore+PDdutycyclezscore+LPstartphasezscore+PYstartphasezscore)/5;
				fitness += 1/(average);
			}
			
		}
		else{
			cout << "possible multi-periodicity" << endl;
			// cout << "LPstartcount = " << LPstartcount << " ,PYstartcount = " << PYstartcount << endl;
			// NO ORDERING POINTS FOR MULTIPERIODIC
		}
	}
	return fitness;

}

int main(int argc, const char* argv[])
{
    ofstream testfile;
    testfile.open("testfile.dat");

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
	CTRNN Agent(3, nullwindowsize, nulllb, nullub, nullbt, nullwt, 16, 16);

    CTRNN Circuit(3, nullwindowsize, nulllb, nullub, nullbt, nullwt, 16, 16);
    char fname[] = "Quinn.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit;

    ifstream HPin;
    HPin.open("HP_unevolved/HPnull.gn");

    Circuit.SetHPPhenotype(HPin,StepSize);

    
    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;
    for (double time = StepSize; time <= RunDuration; time += StepSize) {
        Circuit.EulerStep(StepSize,false,false);
        // testfile << Circuit.outputs(1) << " " << Circuit.outputs(2) << " " << Circuit.outputs(3) << endl;
    }
    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;
    for (double time = StepSize; time <= RunDuration; time += StepSize) {
        Circuit.EulerStep(StepSize,true,false);
        // testfile << Circuit.outputs(1) << " " << Circuit.outputs(2) << " " << Circuit.outputs(3) << endl;
    }

    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;
    cout << PyloricPerformance(Circuit);

    HPin.close();
    testfile.close();

    return 0;
}