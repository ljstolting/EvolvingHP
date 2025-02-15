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

int num_indivs = 26;

ifstream genomesfile;
ofstream outfile;
char genomesfname[] = "../Pyloric CTRNN Genomes/goodgenomesnormalizedfit.dat";
char outfilefname[] = "../Pyloric CTRNN Genomes/testgoodgenomesnormalizedfit.dat";

int main(){
    genomesfile.open(genomesfname);
    outfile.open(outfilefname);
    if (!genomesfile) {
        cerr << "File not found: " << genomesfname << endl;
        exit(EXIT_FAILURE);
    }

    CTRNN Circuit(3);

    TVector<double>phenotype(1,CTRNNphenotypelen);
    double pyl_fitness =  0;
    for (int indiv=1;indiv<=num_indivs;indiv++){

        // for (int i=1;i<=CTRNNphenotypelen;i++){
        //     genomesfile >> phenotype[i];
        // }
        genomesfile >> phenotype; //all at once version
        // read through the fitness value (pyloric fitness)
    
        genomesfile >> pyl_fitness;
        
        phenotype >> Circuit;
        Circuit.RandomizeCircuitOutput(.5,.5);

        // Run the circuit for an initial transient; HP is off and fitness is not evaluated
            for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                Circuit.EulerStep(StepSize,false);
            }

        double pyloricness = 0;

        pyloricness = PyloricPerformance(Circuit); //HP is on during test (set in other file), but null here

        outfile << phenotype << endl;
        outfile << pyloricness << endl << endl;

    }

    return 0;
}