//------------------------------------------------
// Confirm the pyloric fitness of genomes listed
// in a file (such as the generalist optimization
// file) so that HP fitness in each of the
// corresponding subspaces can be scaled by max
// possible (peak) pyloric performance
//------------------------------------------------
#include "CTRNN.h"
#include "pyloric.h"
#include "VectorMatrix.h"

int N = 3;
int CTRNNphenotypelen = (2*N)+(N*N);
double TransientDur = 1000;

int num_indivs = 1;

ifstream genomesfile;
ofstream outfile;
ofstream burstfile;
char genomesfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/80/pyloriccircuit.ns";
char outfilefname[] = "./Specifically Evolved HP mechanisms/Every Circuit/80/pylorictrajectory_test.ns";
char burstfilefname[] = "./Specifically Evolved HP mechanisms/Every Circuit/80/pyloricbursttimes_test.dat";

// double LPbias = 7;
// double PDbias = 1;

int main(){
    genomesfile.open(genomesfname);
    outfile.open(outfilefname);
    burstfile.open(burstfilefname);
    if (!genomesfile) {
        cerr << "File not found: " << genomesfname << endl;
        exit(EXIT_FAILURE);
    }

    CTRNN Circuit(3);
    genomesfile >> Circuit;
    // cout << Circuit.weights << endl;

    // Circuit.SetNeuronBias(1,LPbias);
    // Circuit.SetNeuronBias(3,PDbias);

    // TVector<double>phenotype(1,CTRNNphenotypelen);
    for (int indiv=1;indiv<=num_indivs;indiv++){
        double pyl_fitness =  0;

        // for (int i=1;i<=CTRNNphenotypelen;i++){
        //     genomesfile >> phenotype[i];
        // }
        // genomesfile >> phenotype; //all at once version
        // read through the fitness value (pyloric fitness)
    
        // genomesfile >> pyl_fitness;
        
        // phenotype >> Circuit;
        Circuit.RandomizeCircuitOutput(.5,.5);

        // Run the circuit for an initial transient; HP is off and fitness is not evaluated
        for (double t = StepSize; t <= TransientDur; t += StepSize) {
            Circuit.EulerStep(StepSize,false);
        }

        // TMatrix<double> OutputHist(1,10000,1,3);
        // ifstream outputhistfile;
        // // outputhistfile.open("./Specifically Evolved HP mechanisms/Every Circuit/0/shiftedoutputs.dat");
        // outputhistfile.open("./Specifically Evolved HP mechanisms/Every Circuit/0/outputtrack.dat");
        // outputhistfile >> OutputHist;
        // cout << OutputHist[1][1] << " " << OutputHist[1][2] << " " << OutputHist[1][3] << endl;

        // TVector<double> rhythm_features(1,8);
        // double pyloricness = 0;

        // BurstTimesfromOutputHist(OutputHist, rhythm_features);
        // pyloricness = PyloricFitFromFeatures(rhythm_features);

        pyl_fitness = PyloricPerformance(Circuit,outfile,burstfile);
        double pyl_fitness_2 = PyloricPerformance(Circuit);
        double pyl_fitness_3 = PyloricPerformance(Circuit);

        cout << pyl_fitness << endl;
        cout << pyl_fitness_2 << endl;
        cout << pyl_fitness_3 << endl;

    }
    genomesfile.close();
    outfile.close();
    burstfile.close();

    return 0;
}