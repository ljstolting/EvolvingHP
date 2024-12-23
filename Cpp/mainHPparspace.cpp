//----------------------------------------------------------------------------
// For a given B1,B3 plane (centered on pyloric solution), systematically
// explore the parameter space of HP mechanisms (LBs,ranges,SWs) that succeed.
// I assume (based on my evolutions) that there will be some planes for which
// it is easy and some for which it is hard or impossible. And intermediates for
// which there are two types of mechanisms that work. 
//----------------------------------------------------------------------------

#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

//INDIVIDUALIZED HP MODE
#include <stdio.h>
#include <cstring>
#include <sys/stat.h>
#define Max_Digits 2

// Task params
const double TransientDurationold = 50; //Seconds with HP off
const double PlasticDuration = 10000; //Seconds with HP running
const int N = 3;
const int CTRNNphenotypelen = (2*N)+(N*N);
const int HPphenotypelen = num*4;
const int num_indivs = 1; //how many genomes in the file (100 for Local Run Mode)

// HP parameter space specifications
const double LB1step = .025;
const double LB3step = .025;
const double range = .1; //assume constant range across neurons
const double Btauval = 150; //right in the middle of evol range
const double SWval = 5; //right in middle of evol range (in seconds)

int main(int argc, const char* argv[])
{
    bool range_encoding = true;

    // 2D GRID OF POINTS
    int resolution = 3;
    TVector<double> par_vals(1,resolution);
    par_vals[1] = -10;
    par_vals[2] = 0;
    par_vals[3] = 10;

    int num_pts = pow(resolution,num);
    TMatrix<double> ptlist(1,num_pts,1,num);
    PointGrid(ptlist,par_vals);

    ifstream ifs;

	CTRNN Circuit(3);

    // Local Run Mode
    // char fname[] = "../Pyloric CTRNN Genomes/goodgenomesnormalizedfit.dat";
    // ifs.open(fname);
    // if (!ifs) {
    //     cerr << "File not found: " << fname << endl;
    //     exit(EXIT_FAILURE);
    // }

    // TVector<double>phenotype(1,CTRNNphenotypelen);
    // double pyl_fitness =  0;

    // //Resuming search after a segfault
    // for(int bar=0;bar<10;bar++){
    //     ifs >> phenotype;
    //     ifs >> pyl_fitness;
    // }

    for (int indiv=0;indiv<num_indivs;indiv++){

        //Define pyloric circuit around which to center the slice

        // Local Run Mode
        // ifs >> phenotype;
        // ifs >> pyl_fitness;
        // phenotype >> Circuit;

        // Parallel Supercomputer Mode
        ifs.open("./pyloriccircuit.ns");
        ifs >> Circuit;

        //Define HP parslice output file
        ofstream HPparspacefile;

        // Local Run Mode
        // char indiv_char[Max_Digits + sizeof(char)];
        // std::sprintf(indiv_char, "%d", indiv);
            
        // char outfile[(Max_Digits) + (sizeof(char)) + 51+15];
        // strcpy(outfile, "./Specifically Evolved HP mechanisms/Every Circuit/");
        // strcat(outfile, indiv_char);
        // strcat(outfile, "/HPparslice.dat");

        // Parallel Supercomputer Mode
        char outfile[] = "./HPparslicerangepoint1.dat";

        HPparspacefile.open(outfile);

        //Define HPs based on position in parspace slice
        for (double LB1=0;LB1<=1;LB1+=LB1step){
            for (double LB3=0;LB3<=1;LB3+=LB3step){
                TVector<double> HPphenotype(1,HPphenotypelen);
                int k = 1;
                for (int i=1;i<=num;i++){
                    HPphenotype[k] = Btauval;
                    k++;
                }
                HPphenotype[k] = LB1;
                k++;
                HPphenotype[k] = LB3;
                k++;
                HPphenotype[k] = LB1 + range;
                k++;
                HPphenotype[k] = LB3 + range;
                k++;
                for (int i=1;i<=num;i++){
                    HPphenotype[k] = SWval;
                    k++;
                }
                Circuit.SetHPPhenotype(HPphenotype,StepSize,true);

                //Check grid of initial points to see how many end up pyloric (will just be counting, not keeping fitness)
                int pyloric_count = 0;
                //IN THE BIAS1,BIAS3 PLANE
                for (int i=1;i<=num_pts;i++){
                    Circuit.SetNeuronBias(1,ptlist(i,1));
                    Circuit.SetNeuronBias(3,ptlist(i,2));

                    //Reset Circuit
                    Circuit.RandomizeCircuitOutput(0.5, 0.5);

                    // Run the circuit for an initial transient; HP is off and fitness is not evaluated
                    for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                        Circuit.EulerStep(StepSize,false);
                    }
                    // Apply plasticity for a period of time
                    for (double t = StepSize; t <= PlasticDuration; t += StepSize) {
                        Circuit.EulerStep(StepSize,true);
                    }
                    // Evaluate whether pyloric
                    double pyloricness;
                    pyloricness = PyloricPerformance(Circuit); //HP is left on during test
                    if (pyloricness >= .3){
                        pyloric_count ++;
                    }
                }
                HPparspacefile << pyloric_count << " ";
            }
            HPparspacefile << endl;   
        }
        HPparspacefile.close();
    }
    ifs.close();

    return 0;
}