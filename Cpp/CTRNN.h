// ***********************************************************
// A class for continuous-time recurrent neural networks
//
// RDB
//  8/94 Created
//  12/98 Optimized integration
//  1/08 Added table-based fast sigmoid w/ linear interpolation
// ************************************************************

// Uncomment the following line for table-based fast sigmoid w/ linear interpolation
//#define FAST_SIGMOID

#include "VectorMatrix.h"
#include "random.h"
#include <iostream>
#include <math.h>

#pragma once


// The sigmoid function

#ifdef FAST_SIGMOID
const int SigTabSize = 400;
const double SigTabRange = 15.0;

double fastsigmoid(double x);
#endif

inline double sigma(double x)
{
  return 1/(1 + exp(-x));
}

inline double sigmoid(double x)
{
#ifndef FAST_SIGMOID
  return sigma(x);
#else
  return fastsigmoid(x);
#endif
}


// The inverse sigmoid function

inline double InverseSigmoid(double y)
{
  return log(y/(1-y));
}


// The CTRNN class declaration

class CTRNN {
    public:
        // The constructor
        CTRNN(int size);
        // The destructor
        ~CTRNN();

        // Accessors
        int CircuitSize(void) {return size;};
        void SetCircuitSize(int newsize);
        double NeuronState(int i) {return states[i];};
        double &NeuronStateReference(int i) {return states[i];};
        void SetNeuronState(int i, double value)
            {states[i] = value;outputs[i] = sigmoid(gains[i]*(states[i] + biases[i]));};
        double NeuronOutput(int i) {return outputs[i];};
        double &NeuronOutputReference(int i) {return outputs[i];};
        void SetNeuronOutput(int i, double value)
            {outputs[i] = value; states[i] = InverseSigmoid(value)/gains[i] - biases[i];};
        double NeuronBias(int i) {return biases[i];};
        void SetNeuronBias(int i, double value) {biases[i] = value;};
        double NeuronGain(int i) {return gains[i];};
        void SetNeuronGain(int i, double value) {gains[i] = value;};
        double NeuronTimeConstant(int i) {return taus[i];};
        void SetNeuronTimeConstant(int i, double value) {taus[i] = value;Rtaus[i] = 1/value;};
        double NeuronExternalInput(int i) {return externalinputs[i];};
        double &NeuronExternalInputReference(int i) {return externalinputs[i];};
        void SetNeuronExternalInput(int i, double value) {externalinputs[i] = value;};
        double ConnectionWeight(int from, int to) {return weights[from][to];};
        void SetConnectionWeight(int from, int to, double value) {weights[from][to] = value;};
        // -- NEW
        double NeuronRho(int i) {return rhos[i];};
        void SetNeuronRho(int i, double value) {rhos[i] = value;};
        double PlasticityLB(int i) {return l_boundary[i];};
        void SetPlasticityLB(int i, double value) {l_boundary[i] = value;};
        double PlasticityUB(int i) {return u_boundary[i];};
        void SetPlasticityUB(int i, double value) {u_boundary[i] = value;};
        double NeuronBiasTimeConstant(int i) {return tausBiases[i];};
        void SetNeuronBiasTimeConstant(int i, double value) {tausBiases[i] = value; RtausBiases[i] = 1/value;};
        double ConnectionWeightTimeConstant(int from, int to) {return tausWeights[from][to];};
        void SetConnectionWeightTimeConstant(int from, int to, double value) {tausWeights[from][to] = value; RtausWeights[from][to] = 1/value;};
        int SlidingWindow(int i) {return windowsize[i];};
        // Built in protections against changing step sizes -- entered SW is always time-based
        void SetSlidingWindow(int i, double windsize, double dt) {windowsize[i]=int(1+(windsize/dt));};
        void SetMaxavg(int i, double a) {maxavg[i] = a;};
        void SetMinavg(int i, double a) {minavg[i] = a;};
        // --
        void LesionNeuron(int n)
        {
            for (int i = 1; i<= size; i++) {
                SetConnectionWeight(i,n,0);
                SetConnectionWeight(n,i,0);
            }
        }
        void SetCenterCrossing(void);
        void SetHPPhenotype(istream& is, double dt);
        void SetHPPhenotype(TVector<double>& phenotype, double dt);
        void WriteHPGenome(ostream& os);
        void PrintMaxMinAvgs(void);

        // Input and output
        friend ostream& operator<<(ostream& os, CTRNN& c);
        friend istream& operator>>(istream& is, CTRNN& c);

        // Control
        void RandomizeCircuitState(double lb, double ub);
        void RandomizeCircuitState(double lb, double ub, RandomState &rs);
        void RandomizeCircuitOutput(double lb, double ub);
        void RandomizeCircuitOutput(double lb, double ub, RandomState &rs);
        void RhoCalc(void);
        void EulerStep(double stepsize, bool adaptbiases, bool adaptweights);
        // void RK4Step(double stepsize);

        int size, stepnum;
        TVector<int> windowsize; // NEW for AVERAGING
        double wr, br; // NEWER for CAPPING
        int max_windowsize;
        TVector<double> states, outputs, biases, gains, taus, Rtaus, externalinputs;
        TVector<double> rhos, tausBiases, RtausBiases, l_boundary, u_boundary, minavg, maxavg; // NEW
        TVector<double> avgoutputs; // NEW for AVERAGING
        TMatrix<double> weights;
        TMatrix<double> tausWeights, RtausWeights; // NEW
        TMatrix<double> outputhist; // NEW for AVERAGING
        TVector<double> TempStates,TempOutputs;
};
