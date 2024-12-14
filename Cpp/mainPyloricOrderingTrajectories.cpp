//------------------------------------------------
// For each pyloric genome, collect a multi-cycle
// trajectory for all neurons, and the detected 
// burst starts/ends for future reference
//------------------------------------------------
#include "CTRNN.h"
#include "pyloric.h"
#include "VectorMatrix.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <sys/stat.h> 
#define Max_Digits 10

using namespace std;

int N = 3;
int CTRNNphenotypelen = (2*N)+(N*N);

int num_indivs = 100;

ifstream genomesfile;
ofstream trajfile;
ofstream burstfile;

char genomesfname[] = "../Pyloric CTRNN Genomes/goodgenomesnormalizedfit.dat"; //containing the genomes to be collected and their associated fitnesses
char outdirectory[] = "./Specifically Evolved HP mechanisms/Every Circuit/"; //to be appended to the beginning of each outfile name
char trajfilefname[] = "/pylorictrajectory.ns"; //store the neuron output space trajectory
char burstfilefname[] = "/pyloricbursttimes.dat"; //store the start and end time or each 


int main(){
    genomesfile.open(genomesfname);

    if (!genomesfile) {
        cerr << "File not found: " << genomesfname << endl;
        exit(EXIT_FAILURE);
    }

    CTRNN Circuit(3);

    TVector<double>phenotype(1,CTRNNphenotypelen);
    double pyl_fitness =  0;

    for (int indiv=0;indiv<num_indivs;indiv++){
        char n_char[Max_Digits + sizeof(char)];
        std::sprintf(n_char, "%d", indiv);
            
        char trajfolder[Max_Digits + sizeof(char) + 51+21];
        strcpy(trajfolder, outdirectory);
        strcat(trajfolder, n_char);
        strcat(trajfolder,trajfilefname);
        trajfile.open(trajfolder);

        char burstfolder[Max_Digits + sizeof(char) + 52+22];
        strcpy(burstfolder, outdirectory);
        strcat(burstfolder, n_char);
        strcat(burstfolder,burstfilefname);
        burstfile.open(burstfolder);

        genomesfile >> phenotype; //all at once version
        // read through the fitness value (pyloric fitness)
    
        genomesfile >> pyl_fitness;
        
        phenotype >> Circuit;
        
        // pass transient
        for (double t=StepSize;t<=TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0);
        }

        // record the neurons and the burst properties
        pyl_fitness == PyloricPerformance(Circuit,trajfile,burstfile);

        trajfile.close();
        burstfile.close();
    }

    genomesfile.close();

    return 0;
}