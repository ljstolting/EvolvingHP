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

double TransientDuration = 500;

int N = 3;
int CTRNNphenotypelen = (2*N)+(N*N);

int num_indivs = 100;

ifstream genomesfile;
// ofstream circuitfile;
ifstream circuitfile;
ofstream trajfile;
ofstream burstfile;

// char genomesfname[] = "./pyloric_goodgenomes_timing.dat"; //containing the genomes to be collected and their associated fitnesses
// char outdirectory[] = "./No Timing Requirements/"; //to be appended to the beginning of each outfile name
char outdirectory[] = "./Specifically Evolved HP mechanisms/Every Circuit/";
char circuitfilefname[] = "/pyloriccircuit.ns"; //store (or read) the nervous system
char trajfilefname[] = "/pylorictrajectory_new.dat"; //store the neuron output space trajectory
char burstfilefname[] = "/pyloricbursttimes_new.dat"; //store the start and end time or each 


int main(){
    //FROM FILE INTO SEPARATE FOLDERS MODE
    // genomesfile.open(genomesfname);

    // if (!genomesfile) {
    //     cerr << "File not found: " << genomesfname << endl;
    //     exit(EXIT_FAILURE);
    // }

    CTRNN Circuit(3);

    TVector<double>phenotype(1,CTRNNphenotypelen);
    TVector<double>genotype(1,CTRNNphenotypelen);
    double pyl_fitness =  0;

    for (int indiv=0;indiv<num_indivs;indiv++){
        char n_char[Max_Digits + sizeof(char)];
        std::sprintf(n_char, "%d", indiv);
            
        char trajfolder[Max_Digits + sizeof(char) + 51 + 25];
        strcpy(trajfolder, outdirectory);
        strcat(trajfolder, n_char);
        int result = mkdir(trajfolder,0755);
        strcat(trajfolder,trajfilefname);
        trajfile.open(trajfolder);
        // cout << trajfolder;

        char burstfolder[Max_Digits + sizeof(char) + 51 +26];
        strcpy(burstfolder, outdirectory);
        strcat(burstfolder, n_char);
        strcat(burstfolder,burstfilefname);
        burstfile.open(burstfolder);

        char circuitfolder[Max_Digits + sizeof(char) + 51+18];
        strcpy(circuitfolder, outdirectory);
        strcat(circuitfolder, n_char);
        strcat(circuitfolder,circuitfilefname);
        circuitfile.open(circuitfolder);
        if (!circuitfile) {
            cerr << "File not found: " << circuitfolder << endl;
            exit(EXIT_FAILURE);
        }
        circuitfile.close();
        // cout << circuitfolder << endl;

        //FROM FILE INTO SEPARATE FOLDERS MODE
        // genomesfile >> genotype;
        // genomesfile >> phenotype; //all at once version
        // // read through the fitness value (pyloric fitness)
        // genomesfile >> pyl_fitness;
        
        // phenotype >> Circuit;

        // circuitfile << Circuit;

        //ANALYZE SEPARATE MODE
        cout << circuitfolder << endl;
        circuitfile >> Circuit;
        cout << "testing" << endl;

        Circuit.RandomizeCircuitState(0,0);
        
        // pass transient
        for (double t=StepSize;t<=TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0);
        }

        // record the neurons and the burst properties
        pyl_fitness = PyloricPerformance(Circuit,trajfile,burstfile);

        cout << indiv << ": " << pyl_fitness << endl;;

        trajfile.close();
        burstfile.close();
    }

    // genomesfile.close();

    return 0;
}