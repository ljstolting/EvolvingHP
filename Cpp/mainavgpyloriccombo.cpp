#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

// Task params
const double TransientDuration = 500; //seconds without HP, made pretty long for this analysis (considering HP is off)
const double RunDuration = 50; //seconds to look for an oscillation cycle
const double leaving_tolerance = 0.05; // in state space for greater accuracy
const double return_tolerance = 0.025; //less than leaving tolerance
int max_steps = int((RunDuration*3)/StepSize); //how much memory to allocate for the max size outputhistory

const int N = 3;

// Input files
char resfname[] = "./Test3DHPonPyloricSolutions/res1.dat";  //changes for each compilation
char circuitfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/15/pyloriccircuit.ns";
char dimsfname[] = "./avgsdimensions.dat";

// Output files
char avgsfname[] = "./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice1.dat"; //changes for each compilation
char pylfname[] = "./Test3DHPonPyloricSolutions/pyloricslice3D_15slice1.dat"; //changes for each compilation

int main (int argc, const char* argv[]) 
{
    // Open input files
    ifstream resfile;
	resfile.open(resfname);
	if (!resfile){
		cerr << "File not found: " << resfname << endl;
		exit(EXIT_FAILURE);
	}
	ifstream dimsfile;
	dimsfile.open(dimsfname);
	if (!dimsfile){
		cerr << "File not found: " << dimsfname << endl;
		exit(EXIT_FAILURE);
	}
    ifstream circuitfile;
    circuitfile.open(circuitfname);
    if (!circuitfile) {
        cerr << "File not found: " << circuitfname << endl;
        exit(EXIT_FAILURE);
    }

	// Create files to hold data
	ofstream avgsfile;
	avgsfile.open(avgsfname);
    ofstream pylslicefile;
	pylslicefile.open(pylfname);

    // Organize information about the grid of points, gathered from the res file
    TVector<int> dims(1,N+(N*N));
	int num_dims;

	dimsfile >> dims;
	num_dims = dims.Sum();
	// cout << num_dims;

	TMatrix<double> resmat(1,num_dims,1,3);
	TVector<double> parvec(1,num_dims);

	for(int i=1;i<=num_dims;i++){
		for(int j=1;j<=3;j++){
			resfile >> resmat(i,j);
		}
		parvec(i) = resmat(i,1); //initialize the parameter values at the lowest given bound
	}

    // Create CTRNN and load in parameters
    CTRNN Circuit(N);
    circuitfile >> Circuit;

    Circuit.SetPlasticityPars(dims); //so we can use HP functionality to keep track of the parameters that are being varied

    // Set some variables that will be used in the loop
    TVector<double> acc(1,N); //vector to store the value of the proxy expression for each neuron
    bool finished = false;

    // And finally, start the calculations
	while (!finished){
		// cout << parvec << endl;
        // set the parameters to the new grid value and reset the neural states
		for (int i=1;i<=num_dims;i++){
			Circuit.SetArbDParam(i,parvec(i));
		}
        for (int neuron=1; neuron <= N; neuron ++){
			Circuit.SetNeuronOutput(neuron,.5); //puts in line with the way the average/proxy is calculated most cleanly
		}
		Circuit.WindowReset();

        //pass transient
        for (double t=StepSize;t<=TransientDuration;t+=StepSize){
            Circuit.EulerStep(StepSize,0);
        }
        // then run to keep track of the neuron values either over three full cycles (to ensure three PD starts)
        int cycle_count = 0;
        int stepnum = 0;
        int first_cycle_stepnum = 0;
        TVector<double> startstate(1,N);
        TVector<double> avg(1,N);
        avg.FillContents(0.0);
        TMatrix<double> outputhist(1,max_steps,1,N);
        double t = StepSize; 
        double dist = 0;
        bool left,returned;
        while ((cycle_count<3) && ((t<=RunDuration) || (cycle_count>0)) && (stepnum<max_steps)){
            for (int i = 1; i <= N; i++){
                startstate[i] = Circuit.NeuronState(i);
            }
            // cout << startstate << endl;
            // allow the neural outputs to leave and come back to the same point
            left = false;
            returned = false;
            while ((!returned) && ((t<=RunDuration) || (cycle_count>0)) && (stepnum<max_steps)){
                Circuit.EulerStep(StepSize,0);
                stepnum ++;
                t += StepSize;
                dist = 0;
                for (int i = 1; i <= N; i++){
                    dist += pow(Circuit.NeuronState(i)-startstate[i],2);
                    if (cycle_count == 0){                  //only add up the first cycle
                        avg[i] += Circuit.NeuronOutput(i); 
                    }
                    //record in output hist for pyloric evaluation
                    outputhist(stepnum,i) = Circuit.NeuronOutput(i);
                }
                if (cycle_count == 0){first_cycle_stepnum ++;}

                dist = sqrt(dist);
                // cout << dist << endl;

                if ((dist > leaving_tolerance) && (!left)){
                    left = true;
                }

                if (left && (dist < return_tolerance)){
                    returned = true;
                    left = false;
                    cycle_count ++;
                    // cout << "cycle " << cycle_count << " complete at " << t << endl;
                    // break;
                }
            }
        }
        if (stepnum < max_steps-10){
            //run for 10 more steps just in case we started right on a PD crossing
            for (int s=1;s<=10;s++){
                stepnum ++;
                Circuit.EulerStep(StepSize,0);
                for (int i = 1; i <= N; i++){
                    outputhist(stepnum,i) = Circuit.NeuronOutput(i);
                }
            }
        }
        //take only the part of the output hist that actually ended up getting filled
        TMatrix<double> newoutputhist(1,stepnum,1,N);
        for (int i=1;i<=stepnum;i++){
            for (int j=1;j<=N;j++){
                newoutputhist(i,j) = outputhist(i,j);
            }
        }
        double pyl_fitness;
        if (cycle_count > 0){
            // cout << cycle_count << " cycles" << endl;
            TVector<double> featuresvec(1,8);
            BurstTimesfromOutputHist(newoutputhist,featuresvec);
            pyl_fitness = PyloricFitFromFeatures(featuresvec);
        }
        else{
            pyl_fitness = 0;
        }

        for (int i = 1; i <= N; i++){
            // make into true average by dividing by the number of steps in one cycle
            avg[i] /= first_cycle_stepnum;
        }

        //record all data into files
        avgsfile << avg << endl;
        pylslicefile << pyl_fitness << endl;

        //and then increase the value of the appropriate parameters
		parvec(num_dims)+=resmat(num_dims,3); //step the last dimension
		for (int i=(num_dims-1); i>=1; i-=1){ //start at the second to last dimension and count backwards to see if the next dimension has completed a run
			if(parvec(i+1)>resmat(i+1,2)){   //if the next dimension is over its max
				parvec(i+1) = resmat(i+1,1); //set it to its min
				parvec(i) += resmat(i,3);    //and step the current dimension
                cout << parvec << endl;
			}
		}
		if (parvec(1)>resmat(1,2)){
			finished = true;
		}
    }
    //close relevant files
    resfile.close();
    dimsfile.close();
    circuitfile.close();
    avgsfile.close();
    pylslicefile.close();

    return 0;
}