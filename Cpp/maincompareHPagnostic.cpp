//----------------------------------------------------
// Unfortunately, it appears that there are two anomalies
// that need to be explained. (1) The "HP-agnostic average"
// doesn't seem to correspond with the HP target that would
// produce zero movement at that point. (2) The sliding
// window seems to affect the HP-desired movement, but it
// shouldn't. 
//
// This script is designed to test the space of these anomalies
// with the specific hypothesis that, if you set HP to target
// the HP-agnostic average of a circuit, then it should approach
// no net movement (over time)
//
// It may or may not be the case that it would actually orbit that
// point (timescale sep), but my observations suggest it would
// 
// Note that this is still assuming that statically calculated 
// HP-desired movement  is a good predictor of actual HP movement 
// (timescale separation)
//-----------------------------------------------------
//------------------------------------------------------
// Conclusion drawn: the HP agnostic averages are simply not a
// good predictor of what the HP is going to do at a given point.
// The good news is that the HP desired movement is still a good
// predictor of what HP WILL do at that point. The bad news is I 
// have no idea why that might be.

#include "CTRNN.h"
#include "random.h"
#include "pyloric.h"

using namespace std;

// General run parameters
char CTRNNfname[] = "./Specifically Evolved HP mechanisms/Every Circuit/19/pyloriccircuit.ns";
const double par1min = 13.0595;
const double par1max = 13.0595;
const double par1step = .05;
const double par2min = -15.8667;
const double par2max = -15.8667;
const double par2step = .05;
const double targetstep = .005;

// HP-desired movement parameters (initially taken from other files)
char HPfname[] = "./Convenient HP Mechanisms/bad.dat"; //this file is where you specify the sliding window
                                                           //targets are specified by whatever the agnostic average is
const double RunDuration = 5000;

// HP-agnostic average parameters (initially taken from other files)
const double default_duration = 50;
const double leaving_tolerance = 0.05; // in state space for greater accuracy
const double return_tolerance = 0.025; //less than leaving tolerance

//HP Agnostic Average value at a given point
// Calculated either over one full cycle of the rhythm or the default duration
void HPAgnosticAvg(CTRNN &Agent, TVector<double> &outputvec){
    // Start Circuit from Initial Condition
    Agent.RandomizeCircuitState(.5,.5);

    // Run Transient
    for (double t=StepSize;t<=TransientDuration;t+=StepSize){
        Agent.EulerStep(StepSize,0);
    }

    // then run to keep track of the neuron values either over one full cycle or simply the run duration
    TVector<double> startstate(1,Agent.CircuitSize());
    for (int i = 1; i <= Agent.CircuitSize(); i++){
        startstate[i] = Agent.NeuronState(i);
    }
    // allow the neural outputs to leave and come back to the same point
    bool left = false;
    double dist = 0;
    TVector<double> avg(1,Agent.CircuitSize());
    avg.FillContents(0.0);
    double t = StepSize;
    int stepnum = 0;
    while (t<= RunDuration){
        Agent.EulerStep(StepSize,0);
        stepnum ++;
        t += StepSize;
        dist = 0;
        for (int i = 1; i <= Agent.CircuitSize(); i++){
            // dist += pow(Agent.NeuronState(i)-startstate[i],2); // remove this setting and it will default to always using max time
            avg[i] += Agent.NeuronOutput(i);
        }
        dist = sqrt(dist);

        if ((dist > leaving_tolerance)&&(!left)){
            left = true;
            // cout << "Limit cycle left" << endl;
        }

        if (left && dist < return_tolerance){
            // cout << "Limit cycle came back" << endl;
            break;
        }
    }

    for (int i = 1; i <= Agent.CircuitSize(); i++){
        avg[i] /= stepnum;
        outputvec(i) = avg(i);
    }
    return;
}

void HPDesiredMovement(CTRNN &Agent, TVector<double> &outputvec){
    // Start Circuit from Initial Condition
    Agent.RandomizeCircuitOutput(.5,.5);

    //Run Transient
    for (double t=StepSize;t<=TransientDuration;t+=StepSize){
        Agent.EulerStep(StepSize,0);
    }

    //Keep track of desired bias changes
    for (double t=StepSize;t<=default_duration;t+=StepSize){
        Agent.EulerStepAvgsnoHP(StepSize);
        for (int i = 1; i <= Agent.CircuitSize(); i ++){
            outputvec(i) += StepSize * Agent.RtausBiases(i) * Agent.rhos(i);
        }
    }
    return;
}

int main(int argc, const char* argv[]){
    // double temp = atoi(argv[1]);
    // double bias1fromcall = -10 + temp/20; //will ultimately run for integers between 0 and 600

    // string outfilefname = "./CompareHPagnostic0.dat";

    // outfilefname.replace(19, 1, to_string(atoi(argv[1])));

    string outfilefname = "./CompareHPagnostictest.dat";
    ofstream outfile;
    outfile.open(outfilefname);

    CTRNN Circuit(3);

    ifstream ifs;
    ifs.open(CTRNNfname);
    if (!ifs) {
        cerr << "File not found: " << CTRNNfname << endl;
        exit(EXIT_FAILURE);
    }
    ifs >> Circuit; 
    ifs.close();

    ifstream HPin;
    HPin.open(HPfname);
    if (!HPin) {
        cerr << "File not found: " << HPfname << endl;
        exit(EXIT_FAILURE);
    }
    Circuit.SetHPPhenotype(HPin,StepSize,true);
    HPin.close();

    for (double bias1=par1min;bias1<=par1max;bias1+=par1step){
        cout << bias1 << endl;
        Circuit.SetNeuronBias(1,bias1);
        for (double bias3=par2min;bias3<=par2max;bias3+=par2step){
            Circuit.SetNeuronBias(3,bias3);
            TVector<double> HPagav(1,Circuit.CircuitSize());
            HPAgnosticAvg(Circuit, HPagav);
            cout << "HP Agnostic Averages:" << HPagav << endl;

            // Assuming here that only n1 and n3 are changing
            // uncomment if you are testing difference between given HP and endpoint
            // Circuit.SetPlasticityLB(1,HPagav(1));
            // Circuit.SetPlasticityUB(1,HPagav(1));
            // Circuit.SetPlasticityLB(3,HPagav(3));
            // Circuit.SetPlasticityUB(3,HPagav(3));

            Circuit.SetPlasticityLB(1,.34);
            Circuit.SetPlasticityUB(1,.34);
            Circuit.SetPlasticityLB(3,.34);
            Circuit.SetPlasticityUB(3,.34);
            
            TVector<double> HPdesav(1,Circuit.CircuitSize());
            HPDesiredMovement(Circuit,HPdesav);
            cout << "HP Desired Movement for targets " << Circuit.PlasticityLB(1) << ", " << Circuit.PlasticityLB(3) << ":" << HPdesav << endl;
            outfile << HPdesav(1) << " " << HPdesav(3) << endl;
            // if (HPdesav(1) > 0.01 || HPdesav(3) > 0.01){cout << "large flag" << endl;}

        //     bool N1flip = false;
        //     bool N3flip = false;
        //     TVector<double> HPdesav(1,Circuit.CircuitSize());
        //     HPdesav.FillContents(0);
        //     TVector<double> HPdesav_init(1,Circuit.CircuitSize());
        //     TVector<double> HPdesavhist(1,Circuit.CircuitSize());
        //     TVector<double> flip_targets(1,2);

        //     for (double t=0;t<=1;t+=targetstep){
                
        //         for (int i = 1; i <= 3; i ++){
        //             HPdesavhist(i) = HPdesav(i);
        //         }
                
        //         Circuit.SetPlasticityLB(1,t);
        //         Circuit.SetPlasticityUB(1,t);
        //         Circuit.SetPlasticityLB(3,t); //both neurons can be done at the same time because they're
        //         Circuit.SetPlasticityUB(3,t); // both tracking the same cycle and not changing anything
        //         HPDesiredMovement(Circuit,HPdesav);
        //         // cout << t << " " << HPdesav(1) << " " << HPdesav(3) << " " <<  ((HPdesav(1)*HPdesavhist(1))<0) << " "<<  ((HPdesav(3)*HPdesavhist(3))<0) << endl;
                
        //         if (t == 0){
        //             for (int i=1;i<=Circuit.CircuitSize();i++){
        //                 HPdesav_init(i) = HPdesav(i);
        //             }
        //         }

        //         //record the first positive value of either neurons' desav
        //         if (N1flip == false && (HPdesav(1)*HPdesavhist(1))<0){
        //             flip_targets(1) = t;
        //             N1flip = true;
        //         }

        //         if (N3flip == false && (HPdesav(3)*HPdesavhist(3))<0){
        //             flip_targets(2) = t;
        //             N3flip = true;
        //         }
                
        //         if (N1flip && N3flip){break;}
        //     }
        //     // if one of them did not flip throughout the whole test range, ask whether it was lowest magnitude at 0 or 1
        //     if (!N1flip){
        //         flip_targets(1) = (abs(HPdesav_init(1)) > abs(HPdesav(1))); //using booleans
        //     }
        //     if (!N3flip){
        //         flip_targets(2) = (abs(HPdesav_init(3)) > abs(HPdesav(3)));
            // }
        //     outfile << flip_targets << endl;
        }
        // outfile << endl;
    }
    

    outfile.close();
}