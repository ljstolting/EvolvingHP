//-----------------------------------------------------------
// Run through evolved HP bestind.dat files and collect their
// fitnesses, normalized to be collected in the same way
// despite being from different evoltionary eras
//-----------------------------------------------------------
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"
#include "VectorMatrix.h"
// to do: make the plasticity pars file crutch obsolete

const int Max_Digits = 10;
const bool range_encoding = false;

const int num = 2; //How many parameters does HP have access to change in this data set?

// specify the CTRNN individual on which these HP mechanisms operate
const char Circuitfname[] = "../Pyloric CTRNN Genomes/Pete.ns";
// specify output file name and location
const char outfname[] = "./Specifically Evolved HP mechanisms/Pete/2D/normalizedHPfitsnonrangeencoding.dat";

// STANDARDIZED PARAMETERS THAT AFFECT THE HP FITNESS
// ALREADY DEFINED IN PYLORIC.H (THEREFORE IN THEORY SHOULD BE STANDARDIZED)
// const int num_successive_evaluations = 2; //average the pyloric fitness at how many points in time
// const double TransientDuration = 500; //how long to let the circuit equilibrate before applying HP
const double PlasticDuration = 5000; //how long to allow plasticity to act between the transient and the first evaluation
//                                      //and then between successive evaluations
// const double TestDuration = 100; //maximum number of seconds allowed to test pyloric performance -- can be with HP still on
// const bool HPtest = false; //leave HP on during pyloric evaluation?
// const double StepSize = 0.01 
// const int TestSteps = TestDuration/StepSize; // in steps

const int num_HPs = 212;
const double default_scalingfactor = 25;

double TestFitnessFunction(CTRNN &Agent, TMatrix<double> &ptlist){
	for (int i = 1; i <= ptlist.RowUpperBound();i++){
		Agent.EulerStep(StepSize,true);
	}
	return 15.0;
}

double HPPerformance(CTRNN &Agent, TMatrix<double> &ptlist, double scaling_factor=default_scalingfactor){
	
    double fitness = 0;
	int num_points = ptlist.RowUpperBound();
	// cout << "plane" << Agent.biases << endl;

	for (int i = 1; i <= num_points; i ++){
        // std::cout << "checkpoint 2 " << endl;
		for (int b=1;b<=num;b++){
            // std::cout<< i << " " << b << " " << endl;//<< ptlist(i,b) << endl;
			// Agent.SetArbDParam(b,ptlist(i,b)); 
		}
        Agent.RandomizeCircuitOutput(0.5,0.5);
		// cout << "init " << Agent.biases << endl;
		// cout << "parameters " << Agent.NeuronBias(1) << " " << Agent.NeuronBias(2) << " " << Agent.NeuronBias(3) << " " << Agent.ConnectionWeight(1,1) << endl;
		// Initialize the outputs at 0.5 for all neurons in the circuit
		for (int n=1;n<=Agent.CircuitSize();n++){Agent.SetNeuronState(n,0);}
        // cout << "checkpoint beta " << endl;
		// Run the circuit for an initial transient; HP is off and fitness is not evaluated
		for (double t = StepSize; t <= TransientDuration; t += StepSize) {
			Agent.EulerStep(StepSize,false);
		}
        // cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";

		// Run the circuit for a period of time with HP so the paramters can change
		for (double t = StepSize; t<= PlasticDuration; t+= StepSize){
			Agent.EulerStep(StepSize,true); 
		}
		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << endl;
		double fit = PyloricPerformance(Agent);

		// cout << "fitness" << fit << endl;
		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+scaling_factor;}

		// cout << "scaledfit" << fit;

		// cout << fit << endl;
		fitness += fit;
        // cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << ":";
		// repeat process after another duration of time to ensure solution is stable
		for (double t = StepSize; t<= PlasticDuration; t+= StepSize){
			Agent.EulerStep(StepSize,true);
		}

		// Calculate the Pyloric Fitness
		// cout << Agent.NeuronBias(1) << " " << Agent.NeuronBias(3) << endl;

		fit = PyloricPerformance(Agent);

		// Transform it so Pyloricness at all is worth a lot
		if (fit >= .3){fit = fit+scaling_factor;}

		// cout << "fitness " << i << " " << fit << endl;
		fitness += fit;
		// cout << "final " << Agent.biases << endl;
		// cout << fit << endl << endl;
	}
	// cout << Agent.biases << endl;
    return fitness/(num_points*2);
}

int main(){
// set output file
ofstream outfile;
outfile.open(outfname);

// the initial grid of points from which to evaluate
int resolution = 3;
TVector<double> par_vals(1,resolution);
par_vals[1] = -10;
par_vals[2] = 0;
par_vals[3] = 10;

int num_pts = pow(resolution,num);
TMatrix<double> ptlist(1,num_pts,1,num);
PointGrid(ptlist,par_vals);

// create the CTRNN agent
ifstream Circuitfile;
Circuitfile.open(Circuitfname);
if (!Circuitfile) {
	cerr << "File not found: " << Circuitfname << endl;
	exit(EXIT_FAILURE);
}
CTRNN Agent(3);
Circuitfile >> Agent; 
Circuitfile.close();
 
for (int indiv = 0; indiv < num_HPs; indiv ++){
	// REDEFINING CIRCUIT EVERY TIME AS A CHECK
    // Circuitfile.open(Circuitfname);
	// if (!Circuitfile) {
	// 	cerr << "File not found: " << Circuitfname << endl;
	// 	exit(EXIT_FAILURE);
	// }
	// CTRNN Agent(3);
	// Circuitfile >> Agent; 
	// Circuitfile.close();

	Agent.RandomizeCircuitOutput(.5,.5);

    // import HP

    string bestindfile = "./Specifically Evolved HP mechanisms/Pete/2D/0/bestind.dat";

    bestindfile.replace(45, 1, to_string(indiv));
    // cout << bestindfile << endl;
    ifstream HPbestind;
    HPbestind.open(bestindfile);
    if (!HPbestind) {
        cerr << "File not found: " << bestindfile << endl;
        exit(EXIT_FAILURE);
    }
    Agent.SetHPPhenotype(HPbestind,StepSize,range_encoding);
	// cout << "testing1 " << endl;
	// cout << HPbestind << endl;
	HPbestind.close();
    // cout << "Biases:" << Agent.NeuronBias(1) << endl;
    // cout << "Lower Boundaries:" << Agent.PlasticityLB(1) << " " << Agent.PlasticityLB(3) << endl;
    // cout << Agent.NeuronTimeConstant(1) << endl;
    
    // Evaluate fitness
    double fitness = HPPerformance(Agent,ptlist);
	// double fitness = 10.0;
    // std::cout << "checkpoint charlie" << endl;
    outfile << fitness << endl;
}
// HPbestind.close();

outfile.close();
return 0;
}