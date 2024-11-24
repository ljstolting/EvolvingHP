//------------------------------------------------
// Collect all the pyloric genomes into folders of their own
// to prepare for supercomputer job structure
// Also, collect their normalized fitnesses in a separate file
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
ofstream fitnessfile;
ofstream outfile;

char genomesfname[] = "../Pyloric CTRNN Genomes/goodgenomesnormalizedfit.dat"; //containing the genomes to be collected and their associated fitnesses
char fitnessesfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/normalizedfits.dat"; //collect normalized fitnesses
char outdirectory[] = "./Specifically Evolved HP mechanisms/Every Circuit/"; //to be appended to the beginning of the new directory name
char outfilefname[] = "/pyloriccircuit.ns"; //to be appended to the end of the new directory name

int main(){
    genomesfile.open(genomesfname);
    fitnessfile.open(fitnessesfname);

    if (!genomesfile) {
        cerr << "File not found: " << genomesfname << endl;
        exit(EXIT_FAILURE);
    }
    if (!fitnessfile) {
        cerr << "File not found: " << fitnessesfname << endl;
        exit(EXIT_FAILURE);
    }

    CTRNN Circuit(3);

    TVector<double>phenotype(1,CTRNNphenotypelen);
    double pyl_fitness =  0;

    for (int indiv=0;indiv<num_indivs;indiv++){
        char n_char[Max_Digits + sizeof(char)];
        std::sprintf(n_char, "%d", indiv);
            
        char outfolder[Max_Digits + sizeof(char) + 52+18];
        strcpy(outfolder, outdirectory);
        strcat(outfolder, n_char);
        int result = mkdir(outfolder,0755);

        strcat(outfolder,outfilefname);
        outfile.open(outfolder);

        // for (int i=1;i<=CTRNNphenotypelen;i++){
        //     genomesfile >> phenotype[i];
        // }
        genomesfile >> phenotype; //all at once version
        // read through the fitness value (pyloric fitness)
    
        genomesfile >> pyl_fitness;

        fitnessfile << pyl_fitness << endl; //assume already normalized and calculated
        
        phenotype >> Circuit;
        
        outfile << Circuit;

        outfile.close();
    }

    genomesfile.close();
    fitnessfile.close();

    return 0;
}