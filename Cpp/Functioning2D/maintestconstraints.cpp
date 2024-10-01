// --------------------------------------------------------------
//  Generate random 3N CTRNNs and, if they oscillate, test whether
//  they are guaranteed or unable to satisfy my two derived pyloric
//  conditions. Then, test to see if they are actually pyloric as
//  predicted by the theorem. 
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include <cmath>
#include <charconv>

// Search parameters
const int trials = 100000;  // number or random CTRNNs to generate
const double transient = 500; // seconds
const double max_duration = 100; //seconds given to complete a limit cycle after transient has passed
const double StepSize = .01;
const int max_teststeps = int(max_duration/StepSize);

// Nervous system params
const int N = 3;
const double WR = 16.0; 
const double BR = 16.0; //(WR*N)/2; //<-for allowing center crossing
const double TMIN = .1; 
const double TMAX = 2; 
int	VectSize = N*N + 2*N;
double burstthreshold = 0.5;
double tolerance = .05;

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), -BR, BR);
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				phen(k) = MapSearchParameter(gen(k), -WR, WR);
				k++;
			}
	}
}

int Pyloricness(CTRNN &Agent,int &num_multiperiodic)
{
    string fname;
    fname = "multiperiodicx.dat";

    TMatrix<double> OutputHistory;
	OutputHistory.SetBounds(1,max_teststeps,1,N);
	OutputHistory.FillContents(0.0);
	// TVector<double> CumRateChange(1,N);

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
	while (tstep <= max_teststeps && PDstartcount < 3) {
		// outputhistfile << t << " ";
		for (int i = 1; i <= N; i += 1) {
			OutputHistory(tstep,i) = Agent.NeuronOutput(i);
			// outputhistfile << Agent.NeuronOutput(i) << " ";
			if (Agent.NeuronOutput(i) > maxoutput(i)) {maxoutput(i)=Agent.NeuronOutput(i);}
			if (Agent.NeuronOutput(i) < minoutput(i)) {minoutput(i)=Agent.NeuronOutput(i);}
		}
		// outputhistfile << endl;

		Agent.EulerStep(StepSize,0,0); // noHP

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
				criteriamet += 1;
			}
		}
	}
	if (criteriamet < 3){
		return criteriamet;
	}

	if (PDstartcount < 3){
		cout << "unable to find two full cycles; may want to increase transient, lengthen runtime, or speed up slowest timescale" << endl;
		return criteriamet;
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
		else{cout << "LPstart not found during cycle" << endl; return criteriamet;}

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
		else {cout << "PYstart not found during cycle" << endl; return criteriamet;}


		if ((abs(OutputHistory(PDstarts[1],1) - OutputHistory(PDstarts[2],1))<tolerance)&&(abs(OutputHistory(PDstarts[1],2) - OutputHistory(PDstarts[2],2))<tolerance)){
			// at the two points where PD crosses up, are the other two neurons approximately in the same place?
			// 	ORDERING CRITERIA
			if (LPstart <= PYstart){
				//cout << "order1" << endl;
				criteriamet += 1;
			}
			if (LPend <= PYend){
				//cout << "order2" << endl;
				criteriamet += 1;
			}
			if (PDend <= LPstart){
				//cout << "order3" << endl;
				criteriamet += 1;
			}
        }
		else{
			cout << "possible multi-periodicity " << PDstarts[1] << ": " << OutputHistory(PDstarts[1],1) << ", " << OutputHistory(PDstarts[1],2) << "--" << PDstarts[2] << ": " << OutputHistory(PDstarts[2],1) << ", " << OutputHistory(PDstarts[2],2) << endl;
            // cout << (abs(OutputHistory(PDstarts[1],1) - OutputHistory(PDstarts[2],1))<tolerance) << " " << (abs(OutputHistory(PDstarts[1],2) - OutputHistory(PDstarts[2],2))<tolerance) <<endl;
            // cout << ((abs(OutputHistory(PDstarts[1],1) - OutputHistory(PDstarts[2],1))<tolerance)&&(abs(OutputHistory(PDstarts[1],2) - OutputHistory(PDstarts[2],2))<tolerance)) << endl;
            // // create a file with the output history of this circuit for examination
            // if (num_multiperiodic<10){
            //     ofstream timeseries;
            //     num_multiperiodic ++;
            //     fname[13] = char(48+num_multiperiodic);
            //     timeseries.open(fname);
            //     for (int i=PDstarts[1];i<=PDstarts[3];i++){
            //         timeseries << OutputHistory(i,1) << " " << OutputHistory(i,2) << " " << OutputHistory(i,3) << endl;
            //     }
            //     timeseries.close();
            // }
			// NO ORDERING POINTS FOR MULTIPERIODIC
            return criteriamet;
		}
	}
    
	return criteriamet;
}

bool eq1unsatisfiable(double w12,double w22,double w32,double b2){
    bool cond = min(.5*w12,0.) + min(.5*w32,0.) > -b2 - (.5*w22);
    return cond;
}

bool eq1guaranteed(double w12,double w22,double w32,double b2){
    bool cond = max(.5*w12,0.) + max(.5*w32,0.) < -b2 - (.5*w22);
    return cond;
}

bool eq2unsatisfiable(double w12,double w22,double w32,double b2){
    bool cond = max(.5*w12,w12) + (.5*w22) < min (-b2-w32,-b2);
    return cond;
}

bool eq2guaranteed(double w12,double w22,double w32,double b2){
    bool cond = min(.5*w12,w12) + (.5*w22) > max (-b2-w32,-b2);
    return cond;
}

void RandomizeVector(TVector<double> &v, RandomState rs)
{
	for (int i = 1; i <= v.Size(); i++)
		v[i] = rs.UniformRandom(MinSearchValue,MaxSearchValue); //goes from -1 to 1, hopefully
}

// outfiles
ofstream unsatisfiable; //each row contains status for eq1 and eq2
ofstream guaranteed;    //each row contains status for eq1 and eq2
ofstream pyloric;       //each row contains pyloric status of each already oscillatory circuit
ofstream allcircuits;    //each row contains the number of criteria met for each circuit (3 osc and then 3 ordering if osc)

int main(){
    int num_multiperiodic = 0;
    unsatisfiable.open("unsatisfiable.dat");
    guaranteed.open("guaranteed.dat");
    pyloric.open("pyloric.dat");
    allcircuits.open("allcircuits.dat");
    RandomState rs;
    long IDUM=-time(0);
    int num_osc = 0;
    for (int i=0;i<trials;i++){
        if(i%100==0){cout<< i << " finished, " << num_osc << "oscillating" << endl;}
        // Generate a random genotype
        TVector<double> genotype(1,VectSize);
        IDUM=-time(0);
        rs.SetRandomSeed(IDUM*10*i);
        RandomizeVector(genotype,rs);
        // cout << genotype << endl;

        // Map to phenotype 
        TVector<double> phenotype(1,VectSize);
        GenPhenMapping(genotype,phenotype);

        // Load into nervous system
        TVector<int> Window_Sizes(1,N);
        Window_Sizes.FillContents(1);

        TVector<double> Lower_Bounds(1,N);
        Lower_Bounds.FillContents(0);

        TVector<double> Upper_Bounds(1,N);
        Upper_Bounds.FillContents(1);

        TVector<double> Btaus(1,N);
        Btaus.FillContents(1);

        TMatrix<double> Wtaus(1,N,1,N);
        Wtaus.FillContents(1);

        // Set HP parameters
        CTRNN Circuit(3, Window_Sizes, Lower_Bounds, Upper_Bounds, Btaus, Wtaus, BR, WR);
        int k = 1;
        // Time-constants
        for (int i = 1; i <= N; i++) {
            Circuit.SetNeuronTimeConstant(i,phenotype(k));
            k++;
        }
        // Bias
        for (int i = 1; i <= N; i++) {
            Circuit.SetNeuronBias(i,phenotype(k));
            k++;
        }
        // Weights
        for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    Circuit.SetConnectionWeight(i,j,phenotype(k));
                    k++;
                }
        }
        // Run for transient
        Circuit.RandomizeCircuitOutput(.5,.5); //ignoring multiple attractors
        for (double t=StepSize;t<=transient;t+=StepSize){
            Circuit.EulerStep(StepSize,0,0);
        }
        // Test for (basic) pyloricness
        int pyloric_crit = Pyloricness(Circuit,num_multiperiodic);
        allcircuits << pyloric_crit << endl;
        if (pyloric_crit >= 3){
            num_osc ++;
            bool pyloric_met = (pyloric_crit==6);
            pyloric << pyloric_met << endl;
            // Test parameters for conditions
            unsatisfiable << eq1unsatisfiable(Circuit.weights(1,2),Circuit.weights(2,2),Circuit.weights(3,2),Circuit.biases(2)) << " " << eq2unsatisfiable(Circuit.weights(1,2),Circuit.weights(2,2),Circuit.weights(3,2),Circuit.biases(2)) << endl;
            guaranteed << eq1guaranteed(Circuit.weights(1,2),Circuit.weights(2,2),Circuit.weights(3,2),Circuit.biases(2)) << " " << eq2guaranteed(Circuit.weights(1,2),Circuit.weights(2,2),Circuit.weights(3,2),Circuit.biases(2)) << endl;
        }
    }
    unsatisfiable.close();
    guaranteed.close();
    pyloric.close();
    allcircuits.close();
    return 0;
}