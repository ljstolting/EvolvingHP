import numpy as np
import matplotlib.pyplot as plt

########could someday add option to cut HP off at the designated boundaries or not###############

def sigmoid(x):
    return 1/(1+np.exp(-x))
def invsigmoid(x):
    return np.log(x/(1-x))

#structure of HP genome is:
#[btau1,btau2,btau3,wtau1,wtau2,wtau3,lb1,lb2,lb3,ub1,ub2,ub3,sw1,sw2,sw3]
#assuming each neuron must have same range of acceptability and same temporal averaging of it (SW), but each process (w&b) might have different timescale
class CTRNN():

    def __init__(self,size,dt,duration,HPgenome,neurongenome):
        #define the null HP genome
        if np.all(HPgenome) == None:
            HPgenome = np.ones(5*size)
            HPgenome[size*2:size*3] = 0                          # null lbs = 0
        self.Size = size                                         # number of neurons in the circuit
        self.States = np.ones(size)                              # state of the neurons
        self.Outputs = np.zeros(size)                            # neuron outputs
        self.Weights = np.reshape(neurongenome[0:(size**2)],(size,size))       # weight matrix (numpy array)
        self.Biases = neurongenome[size**2:(size**2)+size]       # bias values (numpy array)
        self.invTimeConstants = 1.0/neurongenome[-size:]         # inverse taus (numpy array)
        self.Inputs = np.zeros((size))                           # external input default to zero
        self.dt = dt                                             # size of integration timestep, in seconds
        self.duration = duration                                 # duration in seconds
        self.time = np.arange(0.0,self.duration,dt)              # timeseries values in seconds
        self.ctrnn_record = np.zeros((size,len(self.time)))      # place to store data of the node outputs over time
        self.Stepnum = 0                                         # initialize the step count at 0

        self.bias_record = np.zeros((size,len(self.time)))       # since parameters of the system are changing under the HP, track biases
        self.weight_record = np.zeros((size,size,len(self.time)))# track weights
        self.invadaptBTimeConsts = 1/HPgenome[0:size]            # bias adaptive time constants 
        self.invadaptWTimeConsts = 1/HPgenome[size:size*2]       # weight adaptive time constants 
        self.lbs = HPgenome[size*2:size*3]                       # lower bounds of the homeostatic target range
        self.ubs = HPgenome[size*3:size*4]                       # upper bounds of the homeostatic target range
        self.slidingwindows = HPgenome[size*4:]                  # how far back to go in timesteps to calculate avg_firingrate (& avg_speed)

        self.rhos = np.zeros(size)                               # store the plastic facilitation parameter of each neuron
        self.max_firingrate = np.zeros(self.Size)                # keep track of the maximum firing rate (for diagnostic)
        self.min_firingrate = np.ones(self.Size)                 # keep track of the minimum firing rate (for diagnostic)
        

    def resetStepcount(self):
        self.Stepnum = 0
        
    def setInputs(self,inputs): #external input to each neuron
        self.Inputs = inputs
    
    def setWeights(self, weights): #weight of connection for each neuron pair, going from row (i) to column (j)
        self.Weights = weights
        
    def randomizeWeights(self):
        self.Weights = np.random.uniform(-16,16,size=(self.Size,self.Size))

    def setBiases(self, biases): #bias shift for each neuron
        self.Biases =  biases
    
    def randomizeBiases(self):
        self.Biases = np.random.uniform(-16,16,size=(self.Size))

    def setTimeConstants(self, timeconstants): #time constant for each neuron
        self.TimeConstants =  np.copy(timeconstants)
        self.invTimeConstants = 1.0/self.TimeConstants
        
    def randomizeTimeConstants(self):
        self.TimeConstants = np.random.uniform(0.5,10,size=(self.Size))
        self.invTimeConstants = 1.0/self.TimeConstants
        
    def setAdaptiveTimeConstants(self, adaptiveWtimeconstants, adaptiveBtimeconstants): #time constants for the adaptation of the weights & biases
        self.invadaptWTimeConsts = 1.0/adaptiveWtimeconstants
        self.invadaptBTimeConsts = 1.0/adaptiveBtimeconstants
        
    def randomizeAdaptiveTimeConstants(self,lower,upper):
        self.invadaptWTimeConsts = 1.0/np.random.uniform(lower,upper,self.Size)
        self.invadaptBTimeConsts = 1.0/np.random.uniform(lower,upper,self.Size)

    def initializeState(self, s):
        self.States = np.copy(s)
        self.Outputs = sigmoid(self.States+self.Biases)

    def initializeOutput(self,o):
        self.Outputs = np.copy(o)
        self.States = invsigmoid(o) - self.Biases
        
    def plasticFacilitationCalc(self): #calculate and update the value of rho for each neuron, using the mean firing rate from the preceding segment of runtime
        for i in range(self.Size):
            if self.Stepnum < self.slidingwindows[i]:
                self.rhos[i] = 0   #not yet enough data to evaluate average firing rate (don't want to use instantaneous firing rate bc it oscillates)
            else:
                avg_firingrate = np.mean(self.ctrnn_record[i,self.Stepnum-int(self.slidingwindows[i]):self.Stepnum])
                if avg_firingrate > self.max_firingrate[i]:
                    self.max_firingrate[i] = avg_firingrate
                if avg_firingrate < self.min_firingrate[i]:
                    self.min_firingrate[i] = avg_firingrate
                if  avg_firingrate < self.lbs[i]:
                    self.rhos[i] = 1-(avg_firingrate/self.lbs[i])
                elif avg_firingrate > self.ubs[i]:
                    self.rhos[i] = (self.ubs[i]-avg_firingrate)/(1-self.ubs[i])
                else:
                    self.rhos[i] = 0  #if in range, no change 
            
    
    def updateBiases(self): #use the value of rho for each neuron to dynamically change biases, scaling the change by 1/speed of walking
        for i in range(self.Size):
            self.Biases[i] += self.dt * self.invadaptBTimeConsts[i] * self.rhos[i]
            if self.Biases[i] >= 16:      #if it goes outside the [-16,16] range, bring it back inside
                self.Biases[i] = 16       #turned off for these experiments
            if self.Biases[i] <= -16:
                self.Biases[i] = -16
    
    def updateWeights(self): #use the value of rho for each neuron to dynamically change all incoming weights to that neuron
        for j in range(self.Size):
            incomingWeights = np.copy(self.Weights[:,j])
            #print("incomingweights to ",j, ": " ,incomingWeights)
            incomingWeights += self.dt * self.invadaptWTimeConsts * self.rhos[j] * np.absolute(incomingWeights)
            self.Weights[:,j] = incomingWeights
            for i in range(self.Size):
                if self.Weights[i,j] >= 16:
                    self.Weights[i,j] = 16
                if self.Weights[i,j] <= -16:
                    self.Weights[i,j] = -16
        
    def ctrnnstep(self,adapt): #use the value of the weights and outputs to change the state of each neuron
        #if adapt = true, then we are implementing the adaptive mechanism
        #first record (initial) condition
        self.Outputs = sigmoid(self.States+self.Biases) 
        self.ctrnn_record[:,self.Stepnum] = self.Outputs #first recording is from the initial condition
        # self.bias_record[:,self.Stepnum]=self.Biases
        # self.weight_record[:,:,self.Stepnum]=self.Weights
        #then update everything
        netinput = self.Inputs + np.dot(self.Weights.T, self.Outputs)
        self.States += self.dt * (self.invTimeConstants*(-self.States+netinput))        
        # self.plasticFacilitationCalc() #moved up here for this notebook
        if adapt == True:
            self.plasticFacilitationCalc()
            self.updateBiases()
            self.updateWeights()
        self.Stepnum += 1
        
    def run(self,adapt):
        for i in range(len(self.time)):
            self.ctrnnstep(adapt)

    def plot(self):
        if self.Size == 3:
            labels = ["LP","PY","PD"]
        else:
            labels = range(self.Size)
        for i in range(self.Size):
            lab = str(labels[i])
            plt.plot(self.time,self.ctrnn_record[i],label=lab)
        plt.plot(self.time,.5*np.ones(len(self.time)))
        plt.title("Neural Activity")
        plt.xlabel("Time (s)")
        plt.ylabel("Firing Rate")
        plt.rcParams["figure.figsize"] = (20,3)
        plt.legend()
        plt.show()
        
    def plotparams(self):
        for i in range(self.Size):
            for j in range(self.Size):
                idx = 3*i+j
                lab = r'$w_{%s%s}$'%(i,j)
                plt.plot(self.time,self.weight_record[i,j,:],label=r"w_{%s%s}"%(i,j))
        for i in range(self.Size):
            plt.plot(self.time,self.bias_record[i,:],label=r'$\theta_%s$'%i)
        plt.plot("CTRNN Parameters")
        plt.xlabel("Time (s)")
        plt.ylabel("Param. Value")
        plt.legend()
        plt.show()

    def pert_step(self,theta,mag):
        '''assumes that the states & outputs have already been intialized on a limit cycle
        theta is in radians'''
        #calculate desired direction of movement
        netinput = self.Inputs + np.dot(self.Weights.T, self.Outputs)
        projected_States = self.States + (self.dt * (self.invTimeConstants*(-self.States+netinput)))
        projected_Outputs = sigmoid(projected_States+self.Biases) 
        projected_direction = projected_Outputs-self.Outputs #in output space
        projected_direction = projected_direction/np.linalg.norm(projected_direction) #normalized, just in case
        #perturb outputs in orthogonal direction, with angle "theta" from the N1 axis
        # first derive rotation matrix https://en.wikipedia.org/wiki/Rotation_matrix
        #then apply it to the reference vector, which RIGHT NOW IS JUST THE Y AXIS PROJECTED UP ONTO THE ORTHOGONAL PLANE
        theta0 = np.array((0,1,projected_direction[1]/projected_direction[2]))
        pert_direction = (np.dot(projected_direction,theta0)*projected_direction)+np.cross(np.cos(theta)*np.cross(projected_direction,theta0),projected_direction)+(np.sin(theta)*np.cross(projected_direction,theta0))
        pert_vector = mag*(pert_direction/np.linalg.norm(pert_direction)) #normalize to 1 and make proper length
        #apply perturbation and align states with perturbation
        self.initializeOutput(self.Outputs+pert_vector)
        #calculate actual direction of movement 
        netinput = self.Inputs + np.dot(self.Weights.T, self.Outputs)
        self.States += self.dt * (self.invTimeConstants*(-self.States+netinput))
        #realign outputs with states and record them
        self.Outputs = sigmoid(self.States+self.Biases) 
        self.ctrnn_record[:,self.Stepnum] = self.Outputs
        #advance simulation
        self.Stepnum += 1