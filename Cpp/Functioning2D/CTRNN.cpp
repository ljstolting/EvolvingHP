// ************************************************************
// HPCTRNN evolution based on Pyloric Fitness
//
// Lindsay Stolting 12/8/22
// ************************************************************

#include "CTRNN.h"
#include "../random.h"
#include "../VectorMatrix.h"
#include <stdlib.h>


// A fast sigmoid implementation using a table w/ linear interpolation
#ifdef FAST_SIGMOID
int SigTableInitFlag = 0;
double SigTab[SigTabSize];

void InitSigmoidTable(void)
{
  if (!SigTableInitFlag) {
    double DeltaX = SigTabRange/(SigTabSize-1);
    for (int i = 0; i <= SigTabSize-1; i++)
      SigTab[i] = sigma(i * DeltaX);
    SigTableInitFlag = 1;
  }
}

double fastsigmoid(double x)
{
  if (x >= SigTabRange) return 1.0;
  if (x < 0) return 1.0 - fastsigmoid(-x);
  double id;
  double frac = modf(x*(SigTabSize-1)/SigTabRange, &id);
  int i = (int)id;
  double y1 = SigTab[i], y2 = SigTab[i+1];

  return y1 + (y2 - y1) * frac;
}
#endif

// ****************************
// Constructors and Destructors
// ****************************

// The constructor
CTRNN::CTRNN(int newsize)
{
	SetCircuitSize(newsize);
#ifdef FAST_SIGMOID
  InitSigmoidTable();
#endif
}

// The destructor

CTRNN::~CTRNN()
{
	SetCircuitSize(0);
}


// *********
// Utilities
// *********

// Resize a circuit.

void CTRNN::SetCircuitSize(int newsize)
{
	size = newsize;
	states.SetBounds(1,size);
	states.FillContents(0.0);
	outputs.SetBounds(1,size);
	outputs.FillContents(0.0);
	biases.SetBounds(1,size);
	biases.FillContents(0.0);
	gains.SetBounds(1,size);
	gains.FillContents(1.0);
	taus.SetBounds(1,size);
	taus.FillContents(1.0);
	Rtaus.SetBounds(1,size);
	Rtaus.FillContents(1.0);
	externalinputs.SetBounds(1,size);
	externalinputs.FillContents(0.0);
	weights.SetBounds(1,size,1,size);
	weights.FillContents(0.0);
	TempStates.SetBounds(1,size);
	TempOutputs.SetBounds(1,size);
  stepnum = 0; 

  // NEW
  rhos.SetBounds(1,size);
  rhos.FillContents(0.0);
  l_boundary.SetBounds(1,size);
  for(int i=1;i<=size;i++){
    l_boundary[i] = 0;
  }
  u_boundary.SetBounds(1,size);
  for(int i=1;i<=size;i++){
    u_boundary[i] = 1;
  }
  tausBiases.SetBounds(1,size);
  RtausBiases.SetBounds(1,size);
  for(int i=1;i<=size;i++){
    tausBiases[i] = 1;
    RtausBiases[i] = 1/tausBiases[i];
  }
  tausWeights.SetBounds(1,size,1,size);
  RtausWeights.SetBounds(1,size,1,size);
  for(int i=1;i<=size;i++){
    for(int j=1;j<=size;j++){
      tausWeights[i][j] = 1;
      RtausWeights[i][j] = 1/(tausWeights[i][j]);
    }
  }

  // NEW for AVERAGING
  windowsize.SetBounds(1,size);
  windowsize.FillContents(1); 
  sumoutputs.SetBounds(1,size);
  sumoutputs.FillContents(0);
  avgoutputs.SetBounds(1,size);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }

  minavg.SetBounds(1,size);
  minavg.FillContents(1);
  maxavg.SetBounds(1,size);
  maxavg.FillContents(0);

  // NEW for CAPPING
  wr = 16;
  br = 16;
}


// *******
// Control
// *******

// Randomize the states or outputs of a circuit.

void CTRNN::RandomizeCircuitState(double lb, double ub)
{
	for (int i = 1; i <= size; i++){
      SetNeuronState(i, UniformRandom(lb, ub));
      SetNeuronOutput(i, sigmoid(gains[i] * (states[i] + biases[i])));
  }
  // reset averaging and sliding window utilities
  minavg.FillContents(1);
  maxavg.FillContents(0);
  sumoutputs.FillContents(0);
  outputhist.FillContents(0.0);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }
  stepnum = 0;
}

void CTRNN::RandomizeCircuitState(double lb, double ub, RandomState &rs)
{
	for (int i = 1; i <= size; i++){
    SetNeuronState(i, rs.UniformRandom(lb, ub));
    SetNeuronOutput(i, sigmoid(gains[i] * (states[i] + biases[i])));
  }
  // reset averaging and sliding window utilities
  minavg.FillContents(1);
  maxavg.FillContents(0);
  sumoutputs.FillContents(0);
  outputhist.FillContents(0.0);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }
  stepnum = 0;
}

void CTRNN::RandomizeCircuitOutput(double lb, double ub)
{
	for (int i = 1; i <= size; i++){
      SetNeuronOutput(i, UniformRandom(lb, ub));
      SetNeuronState(i, (InverseSigmoid(outputs[i])/gains[i])-biases[i]);
  }
  // reset averaging and sliding window utilities
  minavg.FillContents(1);
  maxavg.FillContents(0);
  sumoutputs.FillContents(0);
  outputhist.FillContents(0.0);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }
  stepnum = 0;
}

void CTRNN::RandomizeCircuitOutput(double lb, double ub, RandomState &rs)
{
	for (int i = 1; i <= size; i++){
    SetNeuronOutput(i, rs.UniformRandom(lb, ub));
    SetNeuronState(i, (InverseSigmoid(outputs[i])/gains[i])-biases[i]);
  }
  // reset averaging and sliding window utilities
  minavg.FillContents(1);
  maxavg.FillContents(0);
  sumoutputs.FillContents(0);
  outputhist.FillContents(0.0);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }
  stepnum = 0;
}



// Way to check if all the elements of the output array are now valid CTRNN outputs
// bool checkoutputhist(double array[], int size)
// {
//   for (int i = 0; i < size; i++)
//   {
//       if(array[i] < 0)
//           return false; // return false at the first found

//   }
//   return true; //all elements checked
// }

// Update the averages and rhos of a neuron
void CTRNN::RhoCalc(void){
  // Keep track of the running average of the outputs for some predetermined window of time.
    // 1. Window should always stay updated no matter whether adapting or not (faster so not expensive)
    // 2. Take average for each neuron (unless its sliding window has not yet passed; in that case leave average in between ub and lb to turn HP off)
    for (int i = 1; i <= size; i++){
      // cout << stepnum << " " << windowsize[i] << endl;
      if(stepnum < windowsize[i]){
        outputhist[i][stepnum+1] = NeuronOutput(i);
      }
      if(stepnum == windowsize[i]){ //do initial add-up
      // cout << outputhist(i,1) << endl;
        for (int k = 1; k <= windowsize[i]; k++){  
          sumoutputs[i] += outputhist[i][k];
        }
        // cout << sumoutputs(i) << endl;
        avgoutputs[i] = sumoutputs[i]/windowsize[i];
        if(avgoutputs(i)<minavg(i)){minavg(i)=avgoutputs(i);}; //calc of max and min detected values
        if(avgoutputs(i)>maxavg(i)){maxavg(i)=avgoutputs(i);};
      }
      if(stepnum > windowsize[i]){ //do truncated add-up
        sumoutputs(i) -= outputhist(i,(stepnum%windowsize(i))+1);
        sumoutputs(i) += NeuronOutput(i);
        // replace oldest value
        outputhist(i,(stepnum%windowsize(i))+1) = NeuronOutput(i);
        avgoutputs(i) = sumoutputs(i)/windowsize(i);
        if(avgoutputs(i)<minavg(i)){minavg(i)=avgoutputs(i);}; //calc of max and min detected values
        if(avgoutputs(i)>maxavg(i)){maxavg(i)=avgoutputs(i);};
        // if (i == 2) {cout << stepnum << " " << minavg(2) << " " << maxavg(2) << endl;}
      }
    }


    // NEW: Update rho for each neuron.
    for (int i = 1; i <= size; i++) {
      // cout << l_boundary[i] << " " << u_boundary[i] << endl;
      if (avgoutputs[i] < l_boundary[i]) {
        // cout << l_boundary[i] << endl;
        rhos[i] = (l_boundary[i] - avgoutputs[i])/l_boundary[i];
        // cout << l_boundary[i] << endl << endl;
        // cout << rhos[i] << " " << endl;
      }
      else{
        if (avgoutputs[i] > u_boundary[i]){
          // cout << u_boundary[i] << endl;
          rhos[i] = (u_boundary[i] - avgoutputs[i])/(1.0 - u_boundary[i]);
          // cout << u_boundary[i] << endl << endl;
          // cout << rhos[i] << " " << endl;
        }
        else
        {
          rhos[i] = 0.0; 
        }
      }
      // cout << l_boundary[i] << " " << u_boundary[i] << endl << endl;
    }
}

// Integrate a circuit one step using Euler integration.

void CTRNN::EulerStep(double stepsize, bool adaptbiases, bool adaptweights)
{
  // Update the state of all neurons.
  for (int i = 1; i <= size; i++) {
    double input = externalinputs[i];
    for (int j = 1; j <= size; j++)
      input += weights[j][i] * outputs[j];
    states[i] += stepsize * Rtaus[i] * (input - states[i]);
    outputs[i] = sigmoid(gains[i] * (states[i] + biases[i]));
  }
  if (adaptbiases == true || adaptweights == true){
      RhoCalc();
      stepnum ++;
  }
    // NEW: Update Biases
  if(adaptbiases==true)
  { for (int i = 1; i <= size; i++){
      biases[i] += stepsize * RtausBiases[i] * rhos[i];
      if (biases[i] > br){
          biases[i] = br;
      }
      else{
          if (biases[i] < -br){
              biases[i] = -br;
          }
      }
    } 
  }
  // NEW: Update Weights
  if(adaptweights==true)
  { 
    for (int i = 1; i <= size; i++) 
    {
      for (int j = 1; j <= size; j++)
      {
        weights[i][j] += stepsize * RtausWeights[i][j] * rhos[j] * fabs(weights[i][j]);
        cout << "weight change flag" << endl;
        if (weights[i][j] > wr)
        {
            weights[i][j] = wr;
        }
        else
        {
            if (weights[i][j] < -wr)
            {
                weights[i][j] = -wr;
            }
        }
      }
    }
  }
}

void CTRNN::EulerStepAvgsnoHP(double stepsize)
// Keeps track of the maxmin averages detected, but does not actually change the circuit parameters
{
  // Update the state of all neurons.
  for (int i = 1; i <= size; i++) {
    double input = externalinputs[i];
    for (int j = 1; j <= size; j++)
      input += weights[j][i] * outputs[j];
    states[i] += stepsize * Rtaus[i] * (input - states[i]);
    outputs[i] = sigmoid(gains[i] * (states[i] + biases[i]));
  }
  RhoCalc();
  stepnum ++;
}

void CTRNN::PrintMaxMinAvgs(void){
  cout << "Minimum detected:" << minavg << endl;
  cout << "Maximum detected:" << maxavg << endl;
}



// Set the biases of the CTRNN to their center-crossing values

void CTRNN::SetCenterCrossing(void)
{
    double InputWeights, ThetaStar;

    for (int i = 1; i <= CircuitSize(); i++) {
        // Sum the input weights to this neuron
        InputWeights = 0;
        for (int j = 1; j <= CircuitSize(); j++)
            InputWeights += ConnectionWeight(j, i);
        // Compute the corresponding ThetaStar
        ThetaStar = -InputWeights/2;
        SetNeuronBias(i, ThetaStar);
    }
}

// Define the HP mechanism based on an input file

void CTRNN::SetHPPhenotype(istream& is, double dt, bool range_encoding){
  // Right now, set for the condition where only theta_1 and theta_3 are under HP control
  // Read the bias time constants
  double btau1;
	is >> btau1;
  SetNeuronBiasTimeConstant(1,btau1);

  double btau3;
	is >> btau3;
  SetNeuronBiasTimeConstant(3,btau3);

  // Read the lower bounds
  double lb1;
  is >> lb1;
  SetPlasticityLB(1,lb1);

  double lb3;
  is >> lb3;
  SetPlasticityLB(3,lb3);

  if(range_encoding){
    // Read the ranges and derive the upper bounds
    double range1;
    is >> range1;
    double ub1 = lb1 + range1;
    // if (ub1 > 1){ub1 = 1;} //clipping built into the set funciton
    SetPlasticityUB(1,ub1);

    double range3;
    is >> range3;
    double ub3 = lb3 + range3;
    // if (ub3 > 1){ub3 = 1;} //clipping built into the set funciton
    SetPlasticityUB(3,ub3);
  }

  else
  {  // Read the upper bounds
    double ub1;
    is >> ub1;
    SetPlasticityUB(1,ub1);

    double ub3;
    is >> ub3;
    SetPlasticityUB(3,ub3);
  }

  // Read the sliding windows
  double sw1;
  is >> sw1;
  SetSlidingWindow(1,sw1,dt);

  double sw3;
  is >> sw3;
  SetSlidingWindow(3,sw3,dt);

  // IT IS CRUCIAL TO FIX THE SLIDING WINDOW AVERAGING BEFORE EVALUATION
  // Just in case there is not a transient long enough to fill up the history before HP needs to activate
  max_windowsize = windowsize.Max();
  avgoutputs.SetBounds(1,size);
  for(int i=1;i<=size;i++){
    avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
  }
  outputhist.SetBounds(1,size,1,max_windowsize);
  outputhist.FillContents(0.0);  // fill with zeros 
  // cout << outputhist << endl;
	return;
}

void CTRNN::SetHPPhenotype(TVector<double>& phenotype, double dt, bool range_encoding){
  //Eventually, will work like this
  // int k = 1;
  // for(int i=1;i<=N;i++){
  //     Agent.SetPlasticityBoundary(i,phenotype[k]);
  // }  //....
  SetNeuronBiasTimeConstant(1,phenotype[1]);
  SetNeuronBiasTimeConstant(3,phenotype[2]);
  SetPlasticityLB(1,phenotype[3]);
  SetPlasticityLB(3,phenotype[4]);
  if (range_encoding){
    SetPlasticityUB(1,phenotype[3] + phenotype[5]); 
    SetPlasticityUB(3,phenotype[4] + phenotype[6]);
  }
  else{
  SetPlasticityUB(1,phenotype[5]);
  SetPlasticityUB(3,phenotype[6]);
  }
  SetSlidingWindow(1,phenotype[7],dt); //phenotype sliding window is time based
  SetSlidingWindow(3,phenotype[8],dt);

  // IT IS CRUCIAL TO FIX THE SLIDING WINDOW AVERAGING BEFORE EVALUATION
  max_windowsize = windowsize.Max();
	avgoutputs.SetBounds(1,size);
	for(int i=1;i<=size;i++){
	  avgoutputs[i] = (l_boundary[i]+u_boundary[i])/2; //average of the upper and lower boundaries ensures that initial value keeps HP off
	}
	outputhist.SetBounds(1,size,1,max_windowsize);
	outputhist.FillContents(0.0);  //some number that would never be taken on by the neurons
}

void CTRNN::WriteHPGenome(ostream& os){
  // Right now, set for the condition where only theta_1 and theta_3 are under HP control
  // os << setprecision(32);
  // write the bias time constants
	os << NeuronBiasTimeConstant(1) << " " << NeuronBiasTimeConstant(3) << endl << endl;

  // write the lower bounds
  os << PlasticityLB(1) << " " << PlasticityLB(3) << endl;

  // write the upper bounds
  os << PlasticityUB(1) << " " << PlasticityUB(3) << endl << endl;

  // write the sliding windows
  os << SlidingWindow(1) << " " << SlidingWindow(3);

	return;
}

void CTRNN::SetHPPhenotypebestind(istream &is, double dt, bool range_encoding){
  int trial;
  is >> trial;
  TVector<double> gen(1,4*2); //for now, specific to only two parameters being changed
  TVector<double> phen(1,gen.UpperBound());
  for (int i=1;i<=gen.UpperBound();i++){ 
    is >> gen[i];
  }
  for (int i=1;i<=phen.UpperBound();i++){ 
    is >> phen[i];
  }
  // cout << trial << endl<< gen << endl << phen << endl;
  SetHPPhenotype(phen,dt,range_encoding);
}


// ****************
// Input and Output
// ****************

#include <iomanip>

ostream& operator<<(ostream& os, CTRNN& c)
{//NOT UPDATED TO READ OUT HP PARAMETERS
	// Set the precision
	os << setprecision(32);
	// Write the size
	os << c.size << endl << endl;
	// Write the time constants
	for (int i = 1; i <= c.size; i++)
		os << c.taus[i] << " ";
	os << endl << endl;
	// Write the biases
	for (int i = 1; i <= c.size; i++)
		os << c.biases[i] << " ";
	os << endl << endl;
	// Write the gains
	for (int i = 1; i <= c.size; i++)
		os << c.gains[i] << " ";
	os << endl << endl;
	// Write the weights
	for (int i = 1; i <= c.size; i++) {
		for (int j = 1; j <= c.size; j++)
			os << c.weights[i][j] << " ";
		os << endl;
	}
	// Return the ostream
	return os;
}

istream& operator>>(istream& is, CTRNN& c)
{//NOT UPDATED TO READ IN HP PARAMETERS
	// Read the size
	int size;
	is >> size;
  c.size = size;
	// Read the time constants
	for (int i = 1; i <= size; i++) {
		is >> c.taus[i];
		c.Rtaus[i] = 1/c.taus[i];
	}
	// Read the biases
	for (int i = 1; i <= size; i++)
		is >> c.biases[i];
	// Read the gains
	for (int i = 1; i <= size; i++)
		is >> c.gains[i];
	// Read the weights
	for (int i = 1; i <= size; i++)
		for (int j = 1; j <= size; j++)
			is >> c.weights[i][j];
	// Return the istream
	return is;
}

TVector<double>& operator>>(TVector<double>& phen , CTRNN& c)
{//NOT UPDATED TO READ IN HP PARAMETERS
  int k = 1;
	// Read the time constants
	for (int i = 1; i <= c.size; i++) {
		c.taus[i] = phen[k];
		c.Rtaus[i] = 1/c.taus[i];
    k ++;
	}
	// Read the biases
	for (int i = 1; i <= c.size; i++){
		c.biases[i] = phen[k];
    k ++;
  }
	// Read the weights
	for (int i = 1; i <= c.size; i++){
		for (int j = 1; j <= c.size; j++){
			c.weights[i][j] = phen[k];
      k ++;
    }
  }
	// Return the phenotype
	return phen;
}

TVector<double>& operator<<(TVector<double>& phen, CTRNN& c)
{//NOT UPDATED TO READ OUT HP PARAMETERS -- still CTRNN parameters
  int k = 1;
	// Write the time constants
	for (int i = 1; i <= c.size; i++){
		phen[k] = c.taus[i];
    k++;
  }
	// Write the biases
	for (int i = 1; i <= c.size; i++){
		phen[k] = c.biases[i];
    k ++;
  }

	// Write the weights
	for (int i = 1; i <= c.size; i++) {
		for (int j = 1; j <= c.size; j++){
			phen[k] = c.weights[i][j];
      k ++;
    }
	}
	// Return the phenotype
	return phen;
}
