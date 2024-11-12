//----------------------------------------------------------------------------
// Test a given evolved HP mechanism on all of the evolved pyloric solutions 
// in the "savinggoodgenomes" file
//----------------------------------------------------------------------------

#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

// Task params
const double TransientDurationold = 50; //Seconds with HP off
const double PlasticDuration = 15000; //Seconds with HP running
const int N = 3;
const int CTRNNphenotypelen = (2*N)+(N*N);
const int num_indivs = 90; //how many genomes in the file

int main(int argc, const char* argv[])
{
    ofstream testfile;
    testfile.open("allpyloricXpoint25.dat");

	CTRNN Circuit(3);
    //FILE CONTAINING THE GENOME OF 2-D HP#33, BUT IN THE FORMAT THAT THE CODE READS NOW
    char HPfname[] = "./Functioning2D/HP_unevolved/HPpoint25bestind.gn";
    ifstream HPin;
    HPin.open(HPfname);
    if (!HPin) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    bool range_encoding = true;

    Circuit.SetHPPhenotype(HPin,StepSize,range_encoding);

    // cout << Circuit.PlasticityLB(1) << " " << Circuit.PlasticityLB(2) << " " << Circuit.PlasticityLB(3) << endl;
    // cout << Circuit.plasticitypars << endl << Circuit.plasticneurons << endl; 


    char fname[] = "../Pyloric CTRNN Genomes/goodgenomesnormalizedfit.dat";
    ifstream ifs;
    ifs.open(fname);
    if (!ifs) {
        cerr << "File not found: " << fname << endl;
        exit(EXIT_FAILURE);
    }

    TVector<double>phenotype(1,CTRNNphenotypelen);
    double pyl_fitness =  0;
    int num_perf = 0;
    for (int indiv=1;indiv<=num_indivs;indiv++){

        for (int i=1;i<=CTRNNphenotypelen;i++){
            ifs >> phenotype[i];
        }
        // read through the fitness value (pyloric fitness)
        
        ifs >> pyl_fitness;
        
        phenotype >> Circuit;

        // 2D GRID OF POINTS
        int resolution = 4;
        TVector<double> par_vals(1,resolution);
        par_vals[1] = -16;
        par_vals[2] = -5;
        par_vals[3] = 5;
        par_vals[4] = 16;

        int num_pts = pow(resolution,num);
        TMatrix<double> ptlist(1,num_pts,1,num);
        PointGrid(ptlist,par_vals);

        double pyloricness = 0;
        bool all_pyloric = true; //is true until proven false for one initial condition

        //IN THE BIAS1,BIAS3 PLANE
        for (int i=1;i<=num_pts;i++){
            Circuit.SetNeuronBias(1,ptlist(i,1));
            Circuit.SetNeuronBias(3,ptlist(i,2));
            Circuit.RandomizeCircuitOutput(0.5, 0.5);

            // Run the circuit for an initial transient; HP is off and fitness is not evaluated
            for (double t = StepSize; t <= TransientDuration; t += StepSize) {
                Circuit.EulerStep(StepSize,false);
            }

            for (double t = StepSize; t <= PlasticDuration; t += StepSize) {
                Circuit.EulerStep(StepSize,true);
            }
            pyloricness = PyloricPerformance(Circuit); //HP is left on during test
            testfile << pyloricness << " ";
            if (pyloricness < .3){
                all_pyloric = false;
            }
        }
        num_perf = num_perf + all_pyloric;
        testfile << all_pyloric << endl;
    }

    cout << num_perf;

    HPin.close();
    ifs.close();
    testfile.close();

    return 0;
}