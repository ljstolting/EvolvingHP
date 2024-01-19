#Test the robustness of a pyloric limit cycle to noise
#-------------------------------------------------------------------
#Criteria in paper: 
#LP burst duration, PD burst duration, LP-PD delay, LP-PD interval, PD-LP delay, PD-LP interval, period

#Apply a perturbation of magnitude x in a direciton orthogonal to its direction of movement, then
#calculate resulting limit cycle features. Repeat for m evenly spaced times throughout the original 
#limit cycle, then measure the resulting variability in each criterion by calculating the coefficeint 
#of variation


import numpy as np
from CTRNNclass import *
from pyloricfitness import *

CTRNNsize = 3
HPgenome = np.ones(2*CTRNNsize+3)
HPgenome[0:CTRNNsize] = 0

def pert_run(genome,theta,mag,timelimit,mode = 'bubble',sensitivity_leave=.055,sensitivity_return=.05,plotting=False):
    '''Perturbs in an orthogonal direction at first time step of calculated limit cycle and calculates
    how long it takes to leave and then return within 'sensitivity' to that point in state space. 
    Returns new cycle of the perturbed system restoring itself to the limit cycle and the original limit cycle.'''
    dt=.01
    ctrnn_record = run_timeseries(genome)
    OGperiod, OGlimitcycle = calcfeatures(ctrnn_record)[-2:]
    pert_period = 0
    C = CTRNN(CTRNNsize,dt,timelimit,HPgenome,genome) 
    C.initializeOutput(OGlimitcycle[:,0])
    C.pert_step(theta,mag)
    if mode == 'bubble':
        initial_pt = C.Outputs
        left = 0
        returned = 0
        while (C.Stepnum<=timelimit) and not left:
            C.ctrnnstep(0)
            dist = np.linalg.norm(initial_pt-C.Outputs)
            if dist > sensitivity_leave:
                left = 1
        if left:
            while (C.Stepnum<=timelimit) and not returned:
                C.ctrnnstep(0)
                dist = np.linalg.norm(OGlimitcycle[:,0]-C.Outputs)
                if dist < sensitivity_return:
                    returned = 1
                    pert_period = C.Stepnum
            if not returned:
                print('did not return to bubble')
        else: 
            print('did not leave bubble')
    #mode broken
    # if mode == 'PDstart':
    #     stopped_bursting = 0
    #     cycled=0
    #     former_PDout = 0
    #     while (C.Stepnum<timelimit):
    #         if (C.Outputs[2] < .5) & (former_PDout>= .5):
    #             stopped_bursting=1
    #         if stopped_bursting:
    #             if (C.Outputs[2] > .5) & (former_PDout<= .5):
    #                 pert_period = C.Stepnum
    #                 cycled = 1
    #         former_PDout = C.Outputs[2]
    #         C.ctrnnstep(0)
    #     if not cycled:
    #         print('PD did not initiate burst in given time limit')
    pert_limitcycle = C.ctrnn_record[:,0:pert_period]
    if plotting=='3d':
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot(C.ctrnn_record[0,0:pert_period],C.ctrnn_record[1,0:pert_period],C.ctrnn_record[2,0:pert_period])
        ax.plot(OGlimitcycle[0],OGlimitcycle[1],OGlimitcycle[2],color='orange',alpha=.5)
        plt.show()
        print("new ",pert_period,"old ",OGperiod)
    if plotting=='2d':
        for i in range(3):
            plt.plot(range(OGperiod),OGlimitcycle[i])
        plt.show()
        for i in range(3):
            plt.plot(range(pert_period),pert_limitcycle[i])
        plt.show()
        print("new ",pert_period,"old ",OGperiod)
    return pert_limitcycle, OGlimitcycle

def pert_assessment(genome,sensitivity,timelimit):
    '''Returns average period fluctuation for n perturbations evenly spaced around the limit cycle'''
    return 0

def step_input_run(genome,neuron,n_steps,mag_step,step_len,plotting=False):
    '''go n steps up and then all the way down inputs to specified neuron to output a time series of 
    the perturbed system which can be analyzed'''
    dt=.01 
    pert_mags = np.arange(0,mag_step*n_steps,mag_step)
    pert_mags = np.concatenate((pert_mags,pert_mags[::-1]))
    duration = len(pert_mags)*step_len*dt #in seconds
    C = CTRNN(CTRNNsize,dt,duration,HPgenome,genome)
    C.initializeOutput(np.array([.5,.5,.5]))
    pert_in_form = np.zeros(3)
    pert_in_form[neuron] = 1
    for i in range(n_steps*2):
        C.setInputs(pert_in_form*pert_mags[i])
        for t in range(step_len):
            C.ctrnnstep(0)
    if plotting == True:
        C.plot()
        plt.plot(C.time,np.repeat(pert_mags,step_len))
    return C.ctrnn_record
        
        