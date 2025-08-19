//Update 8/13/25: If the intention of this program is to predict ADHP performance
//  around select bound values, then all data doesn't need to be saved. Only the
//  ongoing successful/mixed/failure status of the points on a more conservative grid. 
//  There is now a new mode that calculates this online instead of collecting all the
//  data and feeding it to a Jupyter notebook. 

#include "TSearch.h"
#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

// Cobbled parallelism
const double slice_step = 0.1;

// Task params
const double TransientDuration = 50; //seconds without HP 
const double RunDuration = 50; //seconds to look for an oscillation cycle
const double leaving_tolerance = 0.05; // in state space for greater accuracy
const double return_tolerance = 0.025; //less than leaving tolerance
int max_steps = int((RunDuration*3)/StepSize); //how much memory to allocate for the max size outputhistory

const int N = 3;

// new mode for keeping track of the files 

// Input files
// char resfname[] = "./Test3DHPonPyloricSolutions/res_high_test.dat";  
// char circuitfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/39/pyloriccircuit.ns";
// char dimsfname[] = "./avgsdimensions.dat";
//Supercomputer
char resfname[] = "../../../Test3DHPonPyloricSolutions/res_high.dat";  
char circuitfname[] = "../../../Specifically Evolved HP mechanisms/Every Circuit/39/pyloriccircuit.ns";
char dimsfname[] = "../../../avgsdimensions.dat";

bool record_all = false;     //previously, we recorded all the averages and their pyloric fitness
                             //for later data analysis (=true), but we can turn this off (=false)
char metaparresfname[] = "../../../metaparres.dat"; //and instead only record the predicted status of the points in HP metapar space
                                             //the resolution of the evaluated points should be the same as the simulated res

// Output files
// LEGACY MODE: track average at each point
// char out1fname[] = "./Test3DHPonPyloricSolutions/HPAgnosticAverage3D_15slice1.dat"; //changes for each compilation
// LEGACY MODE: track pyloric fitness at each point
// char out2fname[] = "./Test3DHPonPyloricSolutions/pyloricslice3D_15slice1.dat"; //changes for each compilation

// NEW MODE: track the predicted HP status at each point in HP metapar space
//  will return two matrices of integers for each point in lattice representing numbers of pyloric
// char out1fname[] = "./Test3DHPonPyloricSolutions/predictedADHPstatus_pyloric_39.dat";
// //  and nonpyloric circuits around that average
// char out2fname[] = "./Test3DHPonPyloricSolutions/predictedADHPstatus_nonpyloric_39.dat"; 

//Supercomputer
char out1fname[] = "./predictedADHPstatus_pyloric_39.dat";
char out2fname[] = "./predictedADHPstatus_nonpyloric_39.dat"; 

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
    ifstream metaparresfile;
    metaparresfile.open(metaparresfname);
    if (!metaparresfile){
        cerr << "File not found:" << metaparresfname << endl;
        exit(EXIT_FAILURE);
    }
    ifstream circuitfile;
    circuitfile.open(circuitfname);
    if (!circuitfile) {
        cerr << "File not found: " << circuitfname << endl;
        exit(EXIT_FAILURE);
    }

	// Create files to hold data
	ofstream outfile1;
	outfile1.open(out1fname);
    ofstream outfile2;
	outfile2.open(out2fname);

    // Organize information about the grid of points, gathered from the res file
    TVector<int> dims(1,N+(N*N));
	int num_dims;

	dimsfile >> dims;
	num_dims = dims.Sum();

	TMatrix<double> resmat(1,num_dims,1,3);
	TVector<double> parvec(1,num_dims);

	for(int i=1;i<=num_dims;i++){
		for(int j=1;j<=3;j++){
			resfile >> resmat(i,j);
		}
	}

    //no longer actually takes first dim from file but takes from the executable call for parallelism
    if (argc>1){
        int slicenum = atoi(argv[1]);
        resmat(1,1) = resmat(1,1)+(slicenum*slice_step); //zero-indexing
        resmat(1,2) = resmat(1,1)+slice_step; //possible that this double counts some stuff, but double counting things won't change the outcome in this scheme
    }   
    cout << resmat(1,1) << " " << resmat(1,2) << endl;

    for(int i=1;i<num_dims;i++){
		parvec(i) = resmat(i,1); //initialize the parameter values at the lowest given bound
    }
    // Create CTRNN and load in parameters
    CTRNN Circuit(N);
    circuitfile >> Circuit;

    Circuit.SetPlasticityPars(dims); //so we can use HP functionality to keep track of the parameters that are being varied
    int regulated_neurons = Circuit.plasticneurons.Sum();
    // cout << "regulated neurons = " << regulated_neurons << endl;

    TMatrix<double> metaparresmat(1,regulated_neurons,1,3);
    for (int i=1;i<=regulated_neurons;i++){
        for (int j=1;j<=3;j++){
            metaparresfile >> metaparresmat(i,j);
            // cout << metaparresmat(i,j) << endl;
        }
    }

    TVector<int> metaparres_lengths(1,regulated_neurons);
    for (int i=1; i<=num_dims; i++){
        int idxlen = (metaparresmat(i,2) - metaparresmat(i,1))/metaparresmat(i,3);
        idxlen = idxlen+1; //upper bound is inclusive
        metaparres_lengths(i) = idxlen;
    }

    //right now, not general to arbitrary dimensions because I don't know enough about the C++ arrays
    int pyloriccount[metaparres_lengths(1)][metaparres_lengths(2)][metaparres_lengths(3)]; 
    int nonpyloriccount[metaparres_lengths(1)][metaparres_lengths(2)][metaparres_lengths(3)];
    // cout << "arrays defined" << endl;
    for (int i=0;i<metaparres_lengths(1);i++){
        for (int j=0;j<metaparres_lengths(2);j++){
            for (int k=0;k<metaparres_lengths(3);k++){
                pyloriccount[i][j][k] = 0;
                nonpyloriccount[i][j][k] = 0;
            }
        }
    }
    // cout << "arrays filled" << endl;

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
        int three_cycle_steps = stepnum;
        while (stepnum < min(three_cycle_steps+10,max_steps)){
            //run for 10 more steps just in case we started right on a PD crossing
            stepnum ++;
            Circuit.EulerStep(StepSize,0);
            for (int i = 1; i <= N; i++){
                outputhist(stepnum,i) = Circuit.NeuronOutput(i);
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

        // LEGACY MODE: record all data into files
        if (record_all){
            outfile1 << avg << endl;
            outfile2 << pyl_fitness << endl;
        }
        // NEW MODE: record the number of pyloric and nonpyloric hp's with average values closest to each point
        else{
            //derive the index of the matrix that corresponds with the point's average values
            TVector<int> add_idx(1,regulated_neurons);
            for (int i=1; i<=regulated_neurons; i++){
                double v = avg(i)-metaparresmat(i,1);
                int idx = round(v/metaparresmat(i,3));
                add_idx(i) = idx;
            }
            if (pyl_fitness>=.3){
                pyloriccount[add_idx(1)][add_idx(2)][add_idx(3)] ++;
            }
            else{
                nonpyloriccount[add_idx(1)][add_idx(2)][add_idx(3)] ++;
            }
        }
        // cout << parvec << endl;
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
            //NEW MODE
            if (!record_all){
                for (int i=0;i<metaparres_lengths(1);i++){
                    for (int j=0;j<metaparres_lengths(2);j++){
                        for (int k=0;k<metaparres_lengths(3);k++){
                            outfile1 << pyloriccount[i][j][k] << " ";
                            outfile2 << nonpyloriccount[i][j][k] << " ";
                        }
                        outfile1 << endl;
                        outfile2 << endl;
                    }
                    outfile1 << endl;
                    outfile2 << endl;
                }
            }
		}
    }
    //close relevant files
    resfile.close();
    dimsfile.close();
    circuitfile.close();
    outfile1.close();
    outfile2.close();

    return 0;
}