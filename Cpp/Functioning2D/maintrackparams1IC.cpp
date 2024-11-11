// --------------------------------------------------------------
//  Track the parameters of a CTRNN as it undergoes Homeostatic Plasticity
// --------------------------------------------------------------
#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"

//#define PRINTOFILE

// Task params
const double TransientDuration = 50; //Seconds with HP off
const double PlasticDuration = 5000; //Seconds with HP running
const double StepSize = 0.01;

//Filenames for parameter traces
char outputswHPfname[] = "../Python/trackoutputscomparewithHP.dat";
char outputswoHPfname[] = "../Python/trackoutputscomparewithoutHP.dat";
char detectedavgsfname[] = "../Python/trackaveragescompare.dat";
char biasesfname[] = "../Python/trackbiasescompare.dat";
// char maxmindebugfname[] = "../Python/maxmindebug.dat";

char HPfname[] = "./18/bestindsfastsuper.dat";

// Will eventually need a function to map genotype of HP into phenotype of values

int main(){
    // Create files to hold data
	ofstream outputswHPfile;
    outputswHPfile.open(outputswHPfname);
    ofstream outputswoHPfile;
    outputswoHPfile.open(outputswoHPfname);
    ofstream avgsfile;
    avgsfile.open(detectedavgsfname);
    ofstream biasesfile;
    biasesfile.open(biasesfname);
    // ofstream maxmindebug;
    // maxmindebug.open(maxmindebugfname);

	// Load the base CTRNN parameters and Set HP parameters
    CTRNN Circuit(3);
    // cout << Circuit.l_boundary << " " << Circuit.u_boundary << endl;
    // cout << Circuit.br;
    char fname[] = "Pete.ns";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
    // Set the proper HP parameters 
    ifstream HPifs;
    HPifs.open(HPfname);
    Circuit.SetHPPhenotypebestind(HPifs,StepSize,false);

    // cout << "Confirm Tau_b's" << Circuit.tausBiases << endl;
    // cout << "Confirm B1: " << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityUB(1) << endl;
    // cout << "Confirm B3: " << Circuit.PlasticityLB(3) << " " << Circuit.PlasticityUB(3) << endl;
    // cout << "Confirm SW: " << Circuit.SlidingWindow(3) << endl;


    // Run a point and record outputs, averages, and biases
    double initialbias1 = 2.59647; //UniformRandom(-16,16);
    double initialbias3 = -8.20154; //UniformRandom(-16,16);
    Circuit.SetNeuronBias(1,initialbias1);
    Circuit.SetNeuronBias(3,initialbias3);
    Circuit.RandomizeCircuitState(0,0);
    for(double t=0;t<TransientDuration;t+=StepSize){
        outputswHPfile << Circuit.outputs << endl;
        biasesfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        // maxmindebug << Circuit.minavg(1) << " " << Circuit.maxavg(1) << " " << Circuit.minavg(3) << " " << Circuit.maxavg(3) << endl;
        avgsfile << Circuit.avgoutputs << endl; //"record" the averages here even though they are just defaults
        Circuit.EulerStep(StepSize,0,0);
    }
    // Circuit.PrintMaxMinAvgs();
    for(double t=0;t<PlasticDuration;t+=StepSize){
        outputswHPfile << Circuit.outputs << endl;
        biasesfile << Circuit.NeuronBias(1) << " " << Circuit.NeuronBias(3) << endl;
        // maxmindebug << Circuit.minavg(1) << " " << Circuit.maxavg(1) << " " << Circuit.minavg(3) << " " << Circuit.maxavg(3) << endl;
        avgsfile << Circuit.avgoutputs << endl;
        Circuit.EulerStep(StepSize,1,0);
    }
    cout << Circuit.biases << endl;

    outputswHPfile.close();
    outputswoHPfile.close();
    biasesfile.close();
    avgsfile.close();

    return 0;
}