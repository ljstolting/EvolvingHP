// ---------------------------------------------------------------
//  Breaking out the family of pyloric fitness functions that are
//  copied in every main file
//
//  Lindsay 4/6/24
// ---------------------------------------------------------------

#include "CTRNN.h"
#include "random.h"
#include "VectorMatrix.h"
#include <iostream>
#include <cmath>

using namespace std;

// Task params
const double TransientDuration = 500; //seconds without HP
const double PlasticDuration1 = 5000; //seconds allowing HP to act
const double PlasticDuration2 = 5000; //number of seconds to wait before testing again, to make sure not relying on precise timing
const double TestDuration = 100; //maximum number of seconds allowed to test pyloric performance -- can be with HP still on
const bool HPtest = true;       //does HP remain on during test (shouldn't matter if platicity time constants are slow enough)
const double StepSize = 0.01;
const int TestSteps = TestDuration/StepSize; // in steps

// Plasticity params 
// MUST MANUALLY CHANGE BECAUSE I CANNOT WORK OUT THE FILE DEPENDENCIES
int num = 3;
int neuronschanging = 3;
int VectSize =  num + (neuronschanging * 3);

// Detection params
const double burstthreshold = .5; //threshold that must be crossed for detecting bursts
const double tolerance = .1; //for detecting double periodicity

double PyloricPerfwTransient(CTRNN &Agent)
{
	int N = Agent.CircuitSize();
	TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);
	double fitness = 0.0;

	// Initialize the outputs at 0.5 for all neurons in the circuit
	Agent.RandomizeCircuitOutput(0.5, 0.5);

	// Run the circuit for an initial transient; HP is off and fitness is not evaluated
	for (double t = StepSize; t <= TransientDuration; t += StepSize) {
		Agent.EulerStep(StepSize,false);
	}

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness while HP is turned OFF.
	int temp = 0;
	for (double time = StepSize; time <= TestDuration; time += StepSize) {
		temp += 1;
		for (int i = 1; i <= N; i += 1) {
			OutputHistory[temp][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,false);
	}
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
	
	if (criteriamet == 3){
		int PDstartcount = 0;
		TVector<int> PDstarts(1,3);
		PDstarts.FillContents(0);
		//LOCATE SECOND TO LAST FULL CYCLE of PD
		for (int step = TestSteps; step >= 2; step --) {
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
			//-----------------
			// cout << LPstart << " " << LPend << " " << PYstart << " " << PYend << " " << PDstarts[1] << " " << PDend << endl;
			//-----------------

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
						// cout << LPstart << ", " << LPend << ", " << PYstart <<", " << PYend << ", " <<PDstarts[1] << ", " <<PDend <<endl;
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

void OrderingRecord(CTRNN &Agent, ofstream &ord_criteria_file){
	//partial copy of PyloricPerf (w/o transient), but prints the fulfilment of ordering criteria to a designated file
	int N = Agent.CircuitSize();
	TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness while HP is turned OFF.
	int temp = 0;
	for (double time = StepSize; time <= TestDuration; time += StepSize) {
		temp += 1;
		for (int i = 1; i <= N; i += 1) {
			OutputHistory[temp][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,false);
	}
	int criteriamet = 0;
	for (int i = 1; i <= N; i += 1) {
		// SHORT HAND FOR ALL NEURONS OSCILLATING APPRECIABLY
		if (minoutput[i] <(burstthreshold-.05)) {
			if (maxoutput[i]>burstthreshold) {
				criteriamet += 1;
			}
		}
	}
	
	if (criteriamet == 3){
		int PDstartcount = 0;
		TVector<int> PDstarts(1,3);
		PDstarts.FillContents(0);
		//LOCATE SECOND TO LAST FULL CYCLE of PD
		for (int step = TestSteps; step >= 2; step --) {
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
			ord_criteria_file << "2 2 2" << endl;
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
						criteriamet += 1;
						ord_criteria_file << 1 << " ";
					}
					else{ord_criteria_file << 0 << " ";}
					if (LPend <= PYend){
						//cout << "order2" << endl;
						criteriamet += 1;
						ord_criteria_file << 1 << " ";
					}
					else{ord_criteria_file << 0 << " ";}
					if (PDend <= LPstart){
						//cout << "order3" << endl;
						criteriamet += 1;
						ord_criteria_file << 1 << endl;
					}
					else{ord_criteria_file << 0 << endl;}
				}
			}
		}
	}
	else{ord_criteria_file << "2 2 2" << endl;}
	// absence of elements printed to file indicates that (1) not all neurons oscillate, (2) not all neurons crossed threshold, or (3) there was "double-periodicity"
	return;
}

double PyloricCountUp(CTRNN &Agent, bool HPtest)
{
	TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,Agent.size);
	OutputHistory.FillContents(0.0);
	double fitness = 0.0;

	TVector<double> maxoutput(1,Agent.size);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,Agent.size);
	minoutput.FillContents(1.0);

	// Run the circuit to calculate Pyloric fitness while HP is on or off depending on input
	int temp = 0;
	for (double time = StepSize; time <= TestDuration; time += StepSize) {
		temp += 1;
		for (int i = 1; i <= Agent.size; i += 1) {
			OutputHistory[temp][i] = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput[i]) {maxoutput[i]=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput[i]) {minoutput[i]=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,HPtest);
	}
	int criteriamet = 0;
	for (int i = 1; i <= Agent.size; i += 1) {
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
		for (int step = TestSteps; step >= 2; step--) {
			if (PDstartcount < 3){
				if (OutputHistory[step][3] > burstthreshold){
					if (OutputHistory[step-1][3] < burstthreshold){
						PDstarts[3-PDstartcount] = step-1;
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
						// cout << LPstart << ", " << LPend << ", " << PYstart <<", " << PYend << ", " <<PDstarts[1] << ", " <<PDend <<endl;
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
        // cout << PDstarts[2] * StepSize << endl;
	}
	return (fitness);

}

double PyloricPerformance(CTRNN &Agent)
//Main difference is that it does not include transient -- assumes circuit is already equilibrated.
//Also only runs until three PDstarts are detected to save time
{
    int N = Agent.CircuitSize();
    TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,TestSteps,1,N);
	OutputHistory.FillContents(0.0);
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
	while (tstep <= TestSteps && PDstartcount < 3) {
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}

		Agent.EulerStep(StepSize,HPtest);

		//Check for PD start
		if (OutputHistory(tstep,3) < burstthreshold && Agent.NeuronOutput(3) > burstthreshold){
			PDstartcount += 1;
			PDstarts[PDstartcount] = tstep;
		}
		tstep += 1;
		t += StepSize;
	}
	for (int i = 1; i <= N; i += 1) {
		// SHORT HAND FOR ALL NEURONS OSCILLATING APPRECIABLY
		if (minoutput[i] <(burstthreshold-.05)) {
			if (maxoutput[i]>burstthreshold) {
				fitness += 0.05;
			}
		}
	}
	if (fitness < 0.15){
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
			}
			if (LPend <= PYend){
				//cout << "order2" << endl;
				fitness += 0.05;
			}
			if (PDend <= LPstart){
				//cout << "order3" << endl;
				fitness += 0.05;
			}
			// cout << LPstart << " " << LPend << " " << PYstart << " " << PYend << " " << PDstarts[1] << " " << PDend << endl;
				
			if (fitness == 0.3){
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
				// cout << LPdutycyclezscore<< ", "<<PYdutycyclezscore<<", "<<PDdutycyclezscore<<", "<<LPstartphasezscore<<", "<<PYstartphasezscore<<endl;
				double average = (LPdutycyclezscore+PYdutycyclezscore+PDdutycyclezscore+LPstartphasezscore+PYstartphasezscore)/5;
				fitness += 1/(average);
			}
			
		}
		else{
			cout << "possible multi-periodicity" << endl;
			// NO ORDERING POINTS FOR MULTIPERIODIC
		}
	}
	return fitness;

}

void converttobase(int N,int resolution,TVector<int> &converted){
	//recursive function to convert to base of choice
	int dimension = converted.UpperBound();
	while(N>0){
		int r = N % resolution;
		N = (N-r) / resolution;
		converted[dimension] = r;
		dimension -= 1;
	}
	return;
}

void PointCombos(TMatrix<int> &answer,int resolution){
	//find all the point combinations for a given dimension count and resolution
	int num_points = answer.RowUpperBound();
	int dimension = answer.ColumnUpperBound();
	for (int i = 1; i <= num_points; i ++){
		TVector<int> row(1,dimension);
		row.FillContents(0);
		converttobase(i-1,resolution,row);
		// cout << row << endl;
		for(int j = 1; j <= dimension; j++){
			answer(i,j) += row(j); //allows for 1 indexing if matrix initialized with ones
		}
	}
	return;
}

double HPPerformance(CTRNN &Agent, double scaling_factor){
	//Starting parameters
	int resolution = 3; //number of points per dimenison

	double lowerres = -8;
	double medres = 0;
	double highres = 8;
	TMatrix<double> par_ICs(1,num,1,resolution);
	// MUST BE CHANGED MANUALLY
	par_ICs.InitializeContents(lowerres,medres,highres,lowerres,medres,highres,lowerres,medres,highres);
	// cout << par_ICs << endl;

	// Agent should already have HP mechanism instantiated
    double fitness = 0;

	int num_points = pow(resolution,num);
	TMatrix<int> par_idxs(1,num_points,1,num);
	par_idxs.FillContents(1);
	
    PointCombos(par_idxs,resolution);
	// cout << par_idxs;

	for (int i = 1; i <= num_points; i ++){
		for (int b=1;b<=num;b++){
			// cout << par_idxs[b] << endl;
			// cout << par_ICs(b,par_idxs(i,b)) << endl;
			Agent.SetNeuronBias(b,par_ICs(b,par_idxs(i,b))); //WILL NEED TO BE GENERALIZED
		}
		// cout << "init" << Agent.biases << endl;
		// cout << "parameters " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
		// Initialize the outputs at 0.5 for all neurons in the circuit
		for (int n=1;n<=Agent.CircuitSize();n++){Agent.SetNeuronState(n,0);}

		// Run the circuit for an initial transient; HP is off and fitness is not evaluated
		for (double t = StepSize; t <= TransientDuration; t += StepSize) {
			Agent.EulerStep(StepSize,false);
		}

		// Run the circuit for a period of time with HP so the paramters can change
		for (double t = StepSize; t<= PlasticDuration1; t+= StepSize){
			Agent.EulerStep(StepSize,true);  //set to only adapt biases, not weights
		}

		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";
		double fit = PyloricPerformance(Agent);

		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+scaling_factor;}

		// cout << fit << endl;
		fitness += fit;

		// repeat process after another duration of time to ensure solution is stable
		for (double t = StepSize; t<= PlasticDuration2; t+= StepSize){
			Agent.EulerStep(StepSize,true);
		}

		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";
		fit = PyloricPerformance(Agent);

		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+scaling_factor;}

		// cout << fit << endl;
		fitness += fit;
		// cout << "final" << Agent.biases << endl;
	
	}
    return fitness/(num_points*2);
}