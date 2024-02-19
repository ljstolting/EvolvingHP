import numpy as np

from CTRNNclass import *

#####To do: make doubly periodic solutions not count##########

#Algorithm parameter adjustments
#----------------------------------------
burst_on_thresh = .5 #the "firing rate threshold" at which the CTRNN neuron is considered to be "bursting"
burst_off_thresh = .45  #the "firing rate" at which the neuron is considered to be effectively silent (just used to determine if oscillation is wide enough)

initial_states = np.array([3.,3.,3.])  #initial states of the neurons
dt=.01
transientdur = 1000 #in seconds
transient = int(transientdur/dt) #in timesteps
plasticdur = 2000 #in seconds
plastic = int(plasticdur/dt)
testdur = 500 #time to simulate CTRNN for in seconds
test = int(testdur/dt)
duration = testdur + plasticdur + transientdur
#----------------------------------------

def run_timeseries(neurongenome,HPgenome = None,plotting=False,debugging=False):
    CTRNNsize = int(np.sqrt(1+len(neurongenome))-1)
    if np.all(HPgenome) == None:
        print('You have not specified an HP genome, so HP will not be used')
        HPgenome = np.ones(2*CTRNNsize+3)
        HPgenome[CTRNNsize*2:CTRNNsize*3] = 0
        HP = 0
    else:
        HP = 1
    C = CTRNN(CTRNNsize,dt,duration,HPgenome,neurongenome)
    C.initializeOutput(np.array([.5,.5,.5]))
    C.resetStepcount()
    for i in range(transient):        #run the CTRNN for the allotted transient duration
        C.ctrnnstep(0)
    for i in range(plastic):        #run the CTRNN for the allotted plasticity duration
        C.ctrnnstep(HP)
    for i in range(test):        #run the CTRNN for the allotted test duration
        C.ctrnnstep(HP)
    if plotting:
        C.plot()
    return C.ctrnn_record

def calcfeatures(CTRNNrecord,debugging=False):
    '''Input is CTRNN time series, uncleansed. 
    Output is start and end for all neurons relative to beginning of cycle, period'''
    #LP = N1, PY = N2, PD = N3
    #Scan to find second to last full PD cycle
    PDstart3 = 0
    PDstart2 = 0
    PDstart1 = 0
    for i in range(len(CTRNNrecord[0]))[:-test:-1]:
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i-1] < burst_on_thresh:
                PDstart3 = i
                break
    for i in range(PDstart3)[:-test:-1]:
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i-1] < burst_on_thresh:
                PDstart2 = i
                break
    for i in range(PDstart2)[:-test:-1]:
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i-1] < burst_on_thresh:
                PDstart1 = i
                break
    if (PDstart1 == 0 or PDstart2 == 0):
        if debugging == True:
            print('unable to find two full cycles,may want to increase runtime')
        return 0,0,0,0,0,0,0,0
    #calculate the start and end times of each neuron's burst in the last full cycle
    LPstart = []
    LPend = 0
    for i in range(PDstart1,PDstart2-1):
        if CTRNNrecord[0,i] < burst_on_thresh:
            if CTRNNrecord[0,i+1] > burst_on_thresh:
                LPstart.append(i)
    if len(LPstart)!=1:
        if debugging == True:
            print('possible double-periodicity')
        return 0,0,0,0,0,0,0,0
    LPstart = LPstart[0]
    for i in range(LPstart,len(CTRNNrecord[0])-1):
        if CTRNNrecord[0,i] > burst_on_thresh:
            if CTRNNrecord[0,i+1] < burst_on_thresh:
                LPend = i
                break
    PYstart = []
    PYend = 0
    for i in range(PDstart1,PDstart2-1):
        if CTRNNrecord[1,i] < burst_on_thresh:
            if CTRNNrecord[1,i+1] > burst_on_thresh:
                PYstart.append(i)
    if len(PYstart)!=1:
        if debugging == True:
            print('possible double-periodicity')
        return 0,0,0,0,0,0,0,0
    PYstart = PYstart[0]
    for i in range(PYstart,len(CTRNNrecord[0])-1):
        if CTRNNrecord[1,i] > burst_on_thresh:
            if CTRNNrecord[1,i+1] < burst_on_thresh:
                PYend = i
                break
    PDend = 0 #end of PD burst
    for i in range(PDstart1,PDstart2): 
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i+1] < burst_on_thresh:
                PDend = i
                break
    limit_cycle = CTRNNrecord[:,PDstart1:PDstart2]
    period=PDstart2-PDstart1
    return LPstart,LPend,PYstart,PYend,PDstart1,PDend,period,limit_cycle

def calcfeatures_accel(CTRNNrecord,debugging=False):
    '''does feature calulation on a pre-cleansed limit cycle, without all the checks for multiperiodicity.
    To be used on perturbed time series, with perturbation at the beginning of the recording. '''
    PDstart1 = 0
    PDstart2 = len(CTRNNrecord[0])-1
    #calculate the start and end times of each neuron's burst in the last full cycle
    LPstart = 0
    LPend = 0
    for i in range(PDstart1,PDstart2-1):
        if CTRNNrecord[0,i] < burst_on_thresh:
            if CTRNNrecord[0,i+1] > burst_on_thresh:
                LPstart = i
                break
    for i in range(LPstart,PDstart2-1):
        if CTRNNrecord[0,i] > burst_on_thresh:
            if CTRNNrecord[0,i+1] < burst_on_thresh:
                LPend = i
                break
    PYstart = 0
    PYend = 0
    for i in range(PDstart1,PDstart2-1):
        if CTRNNrecord[1,i] < burst_on_thresh:
            if CTRNNrecord[1,i+1] > burst_on_thresh:
                PYstart = i
                break
    for i in range(PYstart,PDstart2-1):
        if CTRNNrecord[1,i] > burst_on_thresh:
            if CTRNNrecord[1,i+1] < burst_on_thresh:
                PYend = i
                break
    PDend = 0 #end of PD burst
    for i in range(PDstart1,PDstart2-1): 
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i+1] < burst_on_thresh:
                PDend = i
                break
    period=len(CTRNNrecord[0])
    return LPstart,LPend,PYstart,PYend,PDstart1,PDend,period

def calc_features_subsequent(CTRNNrecord,debugging=False):
    '''give time series with multiple cycles in it and I calculate the features of each cycle'''
    LPstarts = []
    LPends = []
    PYstarts = []
    PYends = []
    PDstarts = []
    PDends = []
    periods = []
    marker = 0
    for i in range(len(CTRNNrecord[0]))[::-1]:
        if CTRNNrecord[2,i] > burst_on_thresh:
            if CTRNNrecord[2,i-1] < burst_on_thresh:
                if marker:
                    PDstarts.append(i)
                    periods.append(formerPDstart-i)
                    formerPDstart = i
                if not marker:
                    formerPDstart = i
                    marker = 1
        if CTRNNrecord[2,i] < burst_on_thresh:
            if CTRNNrecord[2,i-1] > burst_on_thresh:
                if marker:
                    PDends.append(i)
        if CTRNNrecord[1,i] > burst_on_thresh:
            if CTRNNrecord[1,i-1] < burst_on_thresh:
                if marker:
                    PYstarts.append(i)
        if CTRNNrecord[1,i] < burst_on_thresh:
            if CTRNNrecord[1,i-1] > burst_on_thresh:
                if marker:
                    PYends.append(i)
        if CTRNNrecord[0,i] > burst_on_thresh:
            if CTRNNrecord[0,i-1] < burst_on_thresh:
                if marker:
                    LPstarts.append(i)
        if CTRNNrecord[0,i] < burst_on_thresh:
            if CTRNNrecord[0,i-1] > burst_on_thresh:
                if marker:
                    LPends.append(i)
    
    return LPstarts,LPends,PYstarts,PYends,PDstarts,PDends,periods

def pyloriclike_fromfeats(LPstart,LPend,PYstart,PYend,PDstart,PDend,period,debugging=False):
    fitness = .15 #assume all oscillating if using this function
    if LPstart <= PYstart:
            fitness += 0.05
    if LPend <= PYend:
        fitness += 0.05
    if PDend <= LPstart:
        fitness += 0.05
    if fitness==.3:
        LPdutycycle = (LPend-LPstart)/period #burstduration/period
        LPdutycyclezscore = abs(LPdutycycle - .264)/.059
        PYdutycycle = (PYend-PYstart)/period #burstduration/period
        PYdutycyclezscore = abs(PYdutycycle - .348)/.054
        PDdutycycle = (PDend-PDstart)/period #burstduration/period
        PDdutycyclezscore = abs(PDdutycycle - .385)/.040
        LPstartphase = (LPstart-PDstart)/period #delay/period
        LPstartphasezscore = abs(LPstartphase - .533)/.054
        PYstartphase = (PYstart-PDstart)/period #delay/period
        PYstartphasezscore = abs(PYstartphase - .758)/.060
        fitness += 1/(np.average([LPdutycyclezscore,PYdutycyclezscore,PDdutycyclezscore,LPstartphasezscore,PYstartphasezscore]))
        
    if debugging == True:
        print('LPdutycyclezscore ',LPdutycyclezscore)
        print('PYdutycyclezscore ',PYdutycyclezscore)
        print('PDdutycyclezscore ',PDdutycyclezscore)
        print('LPstartphasezscore ',LPstartphasezscore)
        print('PYstartphasezscore ',PYstartphasezscore)
    return fitness

def pyloriclike(neurongenome,HPgenome = None,debugging=False,plotting=False):
    '''input is CTRNN genome [weights,biases,timeconsts] and HP genome is [btau1,btau2,btau3,wtau1,wtau2,wtau3,lb1,lb2,lb3,ub1,ub2,ub3,sw1,sw2,sw3].
    Output is pyloric fitness, with .05 awarded for each oscillating neuron, .05 awarded for each met ordering
    criterion, and additional fitness awarded depending on how close statistics are to experimentally observed
    averages'''
    ctrnn_record = run_timeseries(neurongenome,HPgenome=HPgenome,plotting=plotting,debugging=debugging)
    #check if first three neurons were oscillating (all the way from silent to burst) by the end of the run
    osc = np.zeros(3)
    # print(ctrnn_record[0,-test:])
    for i in range(3):
        #print(max(C.ctrnn_record[i,transient:]),min(C.ctrnn_record[i,transient:]))
        if max(ctrnn_record[i,-test:]) > burst_on_thresh:
            if min(ctrnn_record[i,-test:]) < burst_on_thresh-.025:
                osc[i] = 1
    #print(osc)
    fitness = sum(osc)*0.05 #initialize a fitness value based on how many neurons oscillate sufficiently
    if np.all(osc):
        LPstart,LPend,PYstart,PYend,PDstart,PDend,period = calcfeatures(ctrnn_record[:,-test:])[0:-1]
        if LPstart <= PYstart:
            fitness += 0.05
        if LPend <= PYend:
            fitness += 0.05
        if PDend <= LPstart:
            fitness += 0.05
        if fitness == 0.3:
            LPdutycycle = (LPend-LPstart)/period #burstduration/period
            LPdutycyclezscore = abs(LPdutycycle - .264)/.059
            PYdutycycle = (PYend-PYstart)/period #burstduration/period
            PYdutycyclezscore = abs(PYdutycycle - .348)/.054
            PDdutycycle = (PDend-PDstart)/period #burstduration/period
            PDdutycyclezscore = abs(PDdutycycle - .385)/.040
            LPstartphase = (LPstart-PDstart)/period #delay/period
            LPstartphasezscore = abs(LPstartphase - .533)/.054
            PYstartphase = (PYstart-PDstart)/period #delay/period
            PYstartphasezscore = abs(PYstartphase - .758)/.060
            if debugging == True:
                print('LPdutycyclezscore ',LPdutycyclezscore)
                print('PYdutycyclezscore ',PYdutycyclezscore)
                print('PDdutycyclezscore ',PDdutycyclezscore)
                print('LPstartphasezscore ',LPstartphasezscore)
                print('PYstartphasezscore ',PYstartphasezscore)
            fitness += 1/(np.average([LPdutycyclezscore,PYdutycyclezscore,PDdutycyclezscore,LPstartphasezscore,PYstartphasezscore]))
    return fitness

HPongenome = [.25,.25,.25,.75,.75,.75,40,20,1]
def pyloriclikewithHP(neurongenome,plot):
    return(pyloriclike(neurongenome,HPgenome=HPongenome,plotting=plot))

def continuous_pyloricfitness(neurongenome,HPgenome = None,specificpars=np.ones(15),debugging=False):
    '''New, continuous pyloric fitness function without discrete requirements. Simply does not award extra 
    fitness for oscillation and ordering criteria. Does not check for ordering criteria. If all z-scores were
    met, theoretically ordering criteria would be met, also.'''
    CTRNNsize = int(np.sqrt(1+len(neurongenome))-1)
    if np.all(HPgenome) == None:
        HPgenome = np.ones(2*CTRNNsize+3)
        HPgenome[0:CTRNNsize] = 0
        HP = 0
    else:
        HP = 1
    C = CTRNN(CTRNNsize,dt,duration,HPgenome,neurongenome,specificpars)
    C.initializeState(initial_states)
    C.resetStepcount()
    for i in range(len(C.time)):
        C.ctrnnstep(HP)
    fitness = 0.0
    #check if all neurons were oscillating around the bursting threshold
    osc = np.zeros(3)
    for i in range(3):
        print(max(C.ctrnn_record[i,transient:]),min(C.ctrnn_record[i,transient:]))
        if max(C.ctrnn_record[i,transient:]) > burst_on_thresh+.025:
            if min(C.ctrnn_record[i,transient:]) < burst_on_thresh-.025:
                osc[i] = 1
    if np.all(osc):
        #LP = N1, PY = N2, PD = N3
        #Scan to find second to last full PD cycle
        PDstart3 = 0
        PDstart2 = 0
        PDstart1 = 0
        for i in range(len(C.time))[:transient:-1]:
            if C.ctrnn_record[2,i] > burst_on_thresh:
                if C.ctrnn_record[2,i-1] < burst_on_thresh:
                    PDstart3 = i
                    break
        for i in range(PDstart3)[:transient:-1]:
            if C.ctrnn_record[2,i] > burst_on_thresh:
                if C.ctrnn_record[2,i-1] < burst_on_thresh:
                    PDstart2 = i
                    break
        for i in range(PDstart2)[:transient:-1]:
            if C.ctrnn_record[2,i] > burst_on_thresh:
                if C.ctrnn_record[2,i-1] < burst_on_thresh:
                    PDstart1 = i
                    break
        if (PDstart1 == 0 or PDstart2 == 0):
            if debugging == True:
                print('unable to find two full cycles,may want to increase runtime')
                print('CTRNN',neurongenome)
                #print('HP',HPgenome)
            return fitness
        #calculate the start and end times of each neuron's burst in the last full cycle
        PDend = 0 #end of PD burst
        for i in range(PDstart1,PDstart2): 
            if C.ctrnn_record[2,i] > burst_on_thresh:
                if C.ctrnn_record[2,i+1] < burst_on_thresh:
                    PDend = i
                    break
        LPstart = []
        LPend = 0
        for i in range(PDstart1,PDstart2-1):
            if C.ctrnn_record[0,i] < burst_on_thresh:
                if C.ctrnn_record[0,i+1] > burst_on_thresh:
                    LPstart.append(i)
        if len(LPstart)!=1:
            if debugging == True:
                print('possible double-periodicity')
                print('CTRNN',neurongenome)
                #print('HP',HPgenome)
            return fitness
        LPstart = LPstart[0]
        for i in range(LPstart,len(C.time)-1):
            if C.ctrnn_record[0,i] > burst_on_thresh:
                if C.ctrnn_record[0,i+1] < burst_on_thresh:
                    LPend = i
                    break
        PYstart = []
        PYend = 0
        for i in range(PDstart1,PDstart2-1):
            if C.ctrnn_record[1,i] < burst_on_thresh:
                if C.ctrnn_record[1,i+1] > burst_on_thresh:
                    PYstart.append(i)
        if len(PYstart)!=1:
            if debugging == True:
                print('possible double-periodicity')
                print('CTRNN',neurongenome)
                #print('HP',HPgenome)
            return fitness
        PYstart = PYstart[0]
        for i in range(PYstart,len(C.time)-1):
            if C.ctrnn_record[1,i] > burst_on_thresh:
                if C.ctrnn_record[1,i+1] < burst_on_thresh:
                    PYend = i
                    break
        period = PDstart2 - PDstart1
        LPdutycycle = (LPend-LPstart)/period #burstduration/period
        LPdutycyclezscore = abs(LPdutycycle - .264)/.059
        PYdutycycle = (PYend-PYstart)/period #burstduration/period
        PYdutycyclezscore = abs(PYdutycycle - .348)/.054
        PDdutycycle = (PDend-PDstart1)/period #burstduration/period
        PDdutycyclezscore = abs(PDdutycycle - .385)/.040
        LPstartphase = (LPstart-PDstart1)/period #delay/period
        LPstartphasezscore = abs(LPstartphase - .533)/.054
        PYstartphase = (PYstart-PDstart1)/period #delay/period
        PYstartphasezscore = abs(PYstartphase - .758)/.060
        if debugging == True:
            print('LPdutycyclezscore ',LPdutycyclezscore)
            print('PYdutycyclezscore ',PYdutycyclezscore)
            print('PDdutycyclezscore ',PDdutycyclezscore)
            print('LPstartphasezscore ',LPstartphasezscore)
            print('PYstartphasezscore ',PYstartphasezscore)
        fitness += 1/(np.average([LPdutycyclezscore,PYdutycyclezscore,PDdutycyclezscore,LPstartphasezscore,PYstartphasezscore]))
    return fitness