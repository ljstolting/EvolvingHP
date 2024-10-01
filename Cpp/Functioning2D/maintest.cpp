#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"
#include <chrono>
using namespace std::chrono;

// Task params
const double TransientDurationold = 500; //Seconds with HP off
const double PlasticDuration = 2000; //Seconds with HP running
const double RunDuration = 1000; //How long to test for pyloricness
// const double StepSize = 0.025;
const int RunSteps = RunDuration/StepSize;
const int N = 3;

int main(int argc, const char* argv[])
{
    ofstream testfile;
    testfile.open("testfile.dat");

	CTRNN Circuit(3);
    char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit;

    //incorrectly ordered one
    Circuit.SetNeuronBias(1,10);
    Circuit.SetNeuronBias(3,-5);

	char HPfname[] = "33/bestind.dat";
    ifstream HPin;
    HPin.open(HPfname);
	if (!HPin) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    bool range_encoding = true;

    Circuit.SetHPPhenotypebestind(HPin,StepSize,range_encoding);
    // cout << Circuit.l_boundary(1) << " " << Circuit.l_boundary(3) << endl;
    // cout << Circuit.SlidingWindow(1) << " " << Circuit.SlidingWindow(3) << endl;

	// Initialize states and Run circuit for transient
	Circuit.RandomizeCircuitOutput(0.5, 0.5);

	// Run the circuit for an initial transient; HP is off and fitness is not evaluated
	for (double t = StepSize; t <= TransientDuration; t += StepSize) {
		Circuit.EulerStep(StepSize,false,false);
	}
	auto start = high_resolution_clock::now();
    double scaling_factor = 25;
    
    // for (double t = StepSize; t <= PlasticDuration; t += StepSize) {
	// 	Circuit.EulerStep(StepSize,true,false);
	// }
    cout << "Fitness: " << PyloricPerformance(Circuit) << endl ;
    
    Circuit.PrintMaxMinAvgs();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "duration:" << duration.count() << endl;

	// cout << " " << PyloricPerfwTransient(Circuit);
	// cout << HPPerformance(Circuit,100);

    HPin.close();
    testfile.close();

    return 0;
}