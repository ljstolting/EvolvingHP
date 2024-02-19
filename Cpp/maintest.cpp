#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

const double TransientDuration = 200; //in seconds
const double RunDuration = 200; //in seconds
const double StepSize = 0.01;
const int RunSteps = RunDuration/StepSize; // in steps

const int N = 3;
int	VectSize = N*N + 2*N;

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
    char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit;

    ifstream HPin;
    HPin.open("HP2.gn");

    Circuit.SetHPPhenotype(HPin);

    
    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;
    for (double time = StepSize; time <= RunDuration; time += StepSize) {
        Circuit.EulerStep(StepSize,false,false);
        testfile << Circuit.outputs(1) << " " << Circuit.outputs(2) << " " << Circuit.outputs(3) << endl;
    }
    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;
    for (double time = StepSize; time <= RunDuration; time += StepSize) {
        Circuit.EulerStep(StepSize,true,false);
        testfile << Circuit.outputs(1) << " " << Circuit.outputs(2) << " " << Circuit.outputs(3) << endl;
    }
    // cout << Circuit.NeuronState(1) << endl;
    // cout << Circuit.NeuronBias(1) << endl;

    HPin.close();
    testfile.close();

    return 0;
}