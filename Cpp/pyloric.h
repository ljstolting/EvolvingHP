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
const bool HPtest = false;       //does HP remain on during test (shouldn't matter if platicity time constants are slow enough)
const double StepSize = 0.01;
const int TestSteps = TestDuration/StepSize; // in steps

// Plasticity params 
// MUST MANUALLY CHANGE to reflect plasticpars file BECAUSE I CANNOT WORK OUT THE FILE DEPENDENCIES
int num = 2;
int neuronschanging = 2;
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
// pretty sure this one just records time differently or something
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
	// cout << "checkpoint1" << endl;
    TMatrix<double> OutputHistory(1,TestSteps,1,N);
	// cout << "checkpoint2" << endl;
	
	OutputHistory.FillContents(0.0);

	// cout << "checkpoint3" << endl;
	double fitness = 0.0;

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// cout << "checkpoint4" << endl;

	// Run the circuit to calculate Pyloric fitness -- HP is either left on or turned OFF depending on specification.

	// Run the circuit until you identify 3 PD starts (capping 2 full cycles), keeping track of whether each other neuron crossed the threshold or not
	int tstep = 1;
	double t = StepSize;
	int PDstartcount = 0;
	TVector<int> PDstarts(1,3);
	PDstarts.FillContents(0);
	// cout << "checkpoint5" << endl;

	while (tstep <= TestSteps && PDstartcount < 3) {
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}
		// cout << OutputHistory(1,1) << OutputHistory(1,2) << OutputHistory(1,3);
		// cout << "checkpoint6" << endl;
		Agent.EulerStep(StepSize,HPtest);
		// cout << "checkpoint7" << endl;
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
	// cout << minoutput << maxoutput << endl;
	if (fitness < 0.15){
		// cout << "not oscillatory";
		// cout << "fitness" << fitness;
		return fitness;
	}
	// cout << "checkpoint7" << endl;
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
	// cout << "fitness" << fitness;
	return fitness;

}

//overload but output the neuron trajectories to a file and the burst start/endpoints to another
double PyloricPerformance(CTRNN &Agent, ofstream &trajfile, ofstream &burstfile){
	    int N = Agent.CircuitSize();
	// cout << "checkpoint1" << endl;
    TMatrix<double> OutputHistory(1,TestSteps,1,N);
	// cout << "checkpoint2" << endl;
	
	OutputHistory.FillContents(0.0);

	// cout << "checkpoint3" << endl;
	double fitness = 0.0;

	TVector<double> maxoutput(1,N);
	maxoutput.FillContents(0.0);
	TVector<double> minoutput(1,N);
	minoutput.FillContents(1.0);

	// cout << "checkpoint4" << endl;

	// Run the circuit to calculate Pyloric fitness -- HP is either left on or turned OFF depending on specification.

	// Run the circuit until you identify 3 PD starts (capping 2 full cycles), keeping track of whether each other neuron crossed the threshold or not
	int tstep = 1;
	double t = StepSize;
	int PDstartcount = 0;
	TVector<int> PDstarts(1,3);
	PDstarts.FillContents(0);
	// cout << "checkpoint5" << endl;

	while (tstep <= TestSteps && PDstartcount < 3) {
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}
		Agent.EulerStep(StepSize,HPtest);
		trajfile << Agent.NeuronOutput(1) << " " << Agent.NeuronOutput(2) << " " << Agent.NeuronOutput(3) << endl;

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
	// cout << minoutput << maxoutput << endl;
	if (fitness < 0.15){
		// cout << "not oscillatory";
		// cout << "fitness" << fitness;
		return fitness;
	}
	// cout << "checkpoint7" << endl;
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
				
			if (fitness == 0.3){
				int period = PDstarts[2] - PDstarts[1];

				burstfile << (LPstart-PDstarts[1])*StepSize << " " << (LPend-PDstarts[1])*StepSize << " " << (PYstart-PDstarts[1])*StepSize << " " << (PYend-PDstarts[1])*StepSize << " " << (PDstarts[1]-PDstarts[1])*StepSize << " " << (PDend-PDstarts[1])*StepSize << " " << period*StepSize << endl;

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
	// cout << "fitness" << fitness;
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

void PointGrid(TMatrix<double> &points, TVector<double> &parVals){
	
	//calculate initial conditions on a predefined grid -- gets too unruly when many parameter dimensions are considered
	int resolution = parVals.Size(); //number of points per dimenison

	TMatrix<double> par_ICs(1,num,1,resolution);
	for (int i = 1; i <= num; i++){
		for (int j = 1; j <= resolution; j++){
			par_ICs(i,j) = parVals(j); // assumes that all dimensions have the same resolution of points spaced in the same way
		}
	}

	TMatrix<int>par_idxs(1,points.RowUpperBound(),1,points.ColumnUpperBound());
	par_idxs.FillContents(1);
	
	PointCombos(par_idxs,resolution);
	

	for (int i = 1; i <= points.RowUpperBound(); i ++){
		for (int b=1;b<=num;b++){
			// cout << par_idxs[b] << endl;
			// cout << par_ICs(b,par_idxs(i,b)) << endl;
			points(i,b) = par_ICs(b,par_idxs(i,b));
		}
	}
	return;	
}

