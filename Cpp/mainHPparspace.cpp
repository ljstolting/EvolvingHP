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
const double TransientDuration = 50; //Seconds with HP off
const double PlasticDuration = 100000; //Seconds with HP running
const int N = 3;
const int CTRNNphenotypelen = (2*N)+(N*N);
const int num = 2; //how many paramerters will be under HP control, and perturbed
const int HPphenotypelen = num*4;
const int indiv_idx = 59; //index of the first individual you'll run, only used if file is list
const int num_indivs = 1; //how many genomes in the file (100 for Local Run Mode)
const bool shiftedrho = true;

// HP parameter space specifications
const double LB1min = 0;
const double LB1max = .9;
const double LB3min = 0;
const double LB3max = 0.9;
const double LB1step = .01;
const double LB3step = .01;
const double range = 0.1; //assume constant range across neurons
const double Btauval = 150; //right in the middle of evol range
const double SWval = 0; //(in seconds)

int main(int argc, const char* argv[])
{
    bool range_encoding = true;

    // 2D GRID OF POINTS
    int resolution = 5;
    TVector<double> par_vals(1,resolution);
    par_vals[1] = -10;
    par_vals[2] = -5;
    par_vals[3] = 0;
    par_vals[4] = 5;
    par_vals[5] = 10;

    int num_pts = pow(resolution,num);
    TMatrix<double> ptlist(1,num_pts,1,num);
    PointGrid(ptlist,par_vals);

    //if want from one initial point
    // int num_pts = 1;
    // TMatrix<double> ptlist(1,num_pts,1,num);
    // ptlist(1,1) = -1.67349; 
    // ptlist(1,2) = 1.9336;

    ifstream ifs;

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

    for (int indiv=indiv_idx;indiv<indiv_idx+num_indivs;indiv++){

        //Define pyloric circuit around which to center the slice

        // One circuit only mode
        ifs.open("./Specifically Evolved HP mechanisms/Every Circuit/59/pyloriccircuit.ns");
        CTRNN Circuit(3);
        ifs >> Circuit;

        // Local Run Mode
        // ifs >> phenotype;
        // ifs >> pyl_fitness;
        // phenotype >> Circuit;

        // Parallel Supercomputer Mode
        // ifs.open("./pyloriccircuit.ns");
        // ifs >> Circuit;

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
        // char outfile[] = "./HPparslice_newrho_res5.dat";
        // char outfile[] = "./HPparslice_test.dat";

        // Only one cicuit/HP pair mode
        char outfile[] = "./Specifically Evolved HP mechanisms/Every Circuit/59/HPparslicerangepoint1_rest.dat";

        HPparspacefile.open(outfile);
        ifs.close();

        //Define HPs based on position in parspace slice
        for (double LB1=LB1min;LB1<=LB1max;LB1+=LB1step){
            // cout  << LB1 << endl;
            for (double LB3=LB3min;LB3<=LB3max;LB3+=LB3step){
                // cout << LB1 << " " << LB3 << endl;
                
                Circuit.ShiftedRho(shiftedrho);

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
                HPphenotype[k] = range;
                k++;
                HPphenotype[k] = range;
                k++;
                // cout << HPphenotype << endl;
                for (int i=1;i<=num;i++){
                    HPphenotype[k] = SWval;
                    k++;
                }
                // cout << HPphenotype << Circuit.plasticitypars << endl;
                Circuit.SetHPPhenotype(HPphenotype,StepSize,true);


                cout << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityLB(2) << " " << Circuit.PlasticityLB(3) << " " << Circuit.PlasticityUB(1) << " " << Circuit.PlasticityUB(2) << " " << Circuit.PlasticityUB(3) << endl;
                // cout << "after phen set" << endl;
                //Check grid of initial points to see how many end up pyloric (will just be counting, not keeping fitness)
                int pyloric_count = 0;
                //IN THE BIAS1,BIAS3 PLANE
                for (int i=1;i<=num_pts;i++){
                    // cout << i << endl;
                    Circuit.SetNeuronBias(1,ptlist(i,1));
                    Circuit.SetNeuronBias(3,ptlist(i,2));

                    //Reset Circuit
                    Circuit.RandomizeCircuitOutput(0.5, 0.5);
                    Circuit.WindowReset();

                    // Run the circuit for an initial transient; HP is off and fitness is not evaluated
                    for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                        Circuit.EulerStep(StepSize,false);
                    }
                    
                    // Apply plasticity for a period of time
                    for (double t = StepSize; t <= PlasticDuration; t += StepSize) {
                        Circuit.EulerStep(StepSize,true);
                    }
                    // cout << Circuit.biases << endl;
                    // Evaluate whether pyloric
                    double pyloricness;
                    // cout << "before pyl perf" << endl;
                    pyloricness = PyloricPerformance(Circuit); //HP is off during test?
                    // cout << "after pyl perf" << endl;
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