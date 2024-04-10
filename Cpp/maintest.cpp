#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"
#include <chrono>
using namespace std::chrono;

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

	char HPfname[] = "HP_unevolved/HPtest.gn";
    ifstream HPin;
    HPin.open(HPfname);
	if (!HPin) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }

    Circuit.SetHPPhenotype(HPin,StepSize);

	// Initialize states and Run circuit for transient
	Circuit.RandomizeCircuitOutput(0.5, 0.5);

	// Run the circuit for an initial transient; HP is off and fitness is not evaluated
	for (double t = StepSize; t <= TransientDuration; t += StepSize) {
		Circuit.EulerStep(StepSize,false,false);
	}
	auto start = high_resolution_clock::now();
    cout << PyloricPerformance(Circuit) << endl ;
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