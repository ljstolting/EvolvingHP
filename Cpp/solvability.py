import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

resdata = np.loadtxt('./Specifically Evolved HP mechanisms/Every Circuit/higher_res.dat')
#REGULAR RES
par1min_reg = -16
par1max_reg = 16
par1step_reg = .1
par2min_reg = -16
par2max_reg = 16
par2step_reg = .1

par1min = resdata[0,0]
par1max = resdata[0,1]
par1step = resdata[0,2]
par2min = resdata[1,0]
par2max = resdata[1,1]
par2step = resdata[1,2]

# # which goes an extra for each dimension
par1res = np.arange(par1min,par1max+0.00001,par1step)
par2res = np.arange(par2min,par2max+0.00001,par2step)
par1res_reg = np.arange(par1min_reg,par1max_reg+par1step_reg,par1step_reg)
par2res_reg = np.arange(par2min_reg,par2max_reg+par2step_reg,par2step_reg)
# print(len(par1res))

def get_avg_slice(indiv, type = None):
    if type == "multistability_check":
        data = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/HPAgnosticAverage_highres_multistability.dat"%indiv)
        data = data.reshape((len(par2res),len(par1res),2,3))
        avgs = data[:,:,0,:]
        multistable = np.zeros_like(avgs[:,:,0])
        multistable[np.where(data[:,:,1,0]!=0)] = 1
        avgs = np.swapaxes(avgs,1,0)
        multistable = multistable.swapaxes(0,1)
        return avgs, multistable
    if type == "high_res":
        avgs = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/HPAgnosticAverage_highres.dat"%indiv)
        avgs = avgs.reshape((len(par2res),len(par1res),3))
        avgs = np.swapaxes(avgs,1,0)
        return avgs
    if type == 'netchange':
        avgs = np.loadtxt("./Convenient HP Mechanisms/goodnetchange.dat")
        avgs = avgs.reshape((len(par2res),len(par1res),2))
        avgs = np.swapaxes(avgs,1,0)
        return avgs
    if type == 'newrho':
        data = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/HPAgnosticAverage_highres_newrho.dat"%indiv)
        data = data.reshape((len(par2res),len(par1res),2,3))
        avgs = data[:,:,0,:]
        avgs = np.swapaxes(avgs,1,0)
        dutycycles = data[:,:,1,:]
        dutycycles = dutycycles.swapaxes(1,0)
        return avgs,dutycycles
    else:
        avgs = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/HPAgnosticAverage.dat"%indiv)
        avgs = avgs.reshape((len(par2res_reg),len(par1res_reg),3))
        avgs = np.swapaxes(avgs,1,0)
        return avgs

def get_proxy_slice(indiv):
    '''returns a masked array with the averages of all the stable points, and a list of lambda functions for each oscillatory point, which takes the variables lb and ub and determines HP movement '''
    data = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/HPAgnosticAverage_highres_newrho.dat"%indiv)
    data = data.reshape(len(par2res),len(par1res),2,3)
    avgs = data[:,:,0,:]
    avgs = avgs.swapaxes(0,1)
    duty_cycles = data[:,:,1,:]
    duty_cycles = duty_cycles.swapaxes(0,1)
    masked_stable_avgs = np.ma.masked_where(duty_cycles != 0,avgs)
    masked_oscillatory_avgs = np.ma.masked_where(duty_cycles == 0, avgs)
    masked_oscillatory_dutycycles = np.ma.masked_where(duty_cycles == 0, duty_cycles)

    return masked_stable_avgs, masked_oscillatory_avgs, masked_oscillatory_dutycycles

def get_osc_proxy(masked_oscillatory_avgs, masked_oscillatory_dutycycles,LB,UB):
    pred_HP_movement = (np.ones_like(masked_oscillatory_dutycycles)-masked_oscillatory_dutycycles)*LB + (masked_oscillatory_dutycycles)*UB - masked_oscillatory_avgs
    # pred_HP_movement = gaussian_filter(pred_HP_movement, sigma=.2) #doesn't really work, just hugs edges
    return pred_HP_movement
    

def get_pyloric_slice(indiv,high_res=False):
    if high_res:
        pyloricslice = np.loadtxt('./Specifically Evolved HP mechanisms/Every Circuit/%s/pyloricslice_highres.dat'%indiv).reshape((len(par1res),len(par2res)))
    else:
        pyloricslice = np.loadtxt('./Specifically Evolved HP mechanisms/Every Circuit/%s/pyloricslice.dat'%indiv).reshape((len(par1res_reg),len(par2res_reg)))
    pyloricslice = np.swapaxes(pyloricslice,1,0)
    return pyloricslice

def get_pyloric_outline(indiv,high_res=False):
    pyloricslice = get_pyloric_slice(indiv,high_res)
    borderlist_left = []
    borderlist_right = []
    for i in range(len(pyloricslice)):
        for j in range(2,len(pyloricslice[0])-1):
            if (pyloricslice[i,j-2] < .3 and pyloricslice[i,j-1] < .3 and pyloricslice[i,j]>=.3 and pyloricslice[i,j+1]>=.3):
                borderlist_left.append([par1res[j-1],par2res[i]])
                # print('left @ (%s,%s):'%(par2res[j],par1res[i]),pyloricslice[i,j-1],' right:',pyloricslice[i,j])
            if (pyloricslice[i,j-2] >= .3 and pyloricslice[i,j-1] >= .3 and pyloricslice[i,j]<.3 and pyloricslice[i,j+1]<.3):
                borderlist_right.append([par1res[j],par2res[i]])
    return [borderlist_left,borderlist_right[::-1]]

def uniquevals(indiv,tolerance = 0.03,evaluate = False,plot=True,ax=None,progress=False,high_res=True):
    if high_res:
        stype = 'high_res'
    else:
        stype = None
    pyloricslice = get_pyloric_slice(indiv, high_res)
    avgs = get_avg_slice(indiv,type=stype)
    # maskingslice = np.repeat(pyloricslice,3,axis=1).reshape((len(pyloricslice),len(pyloricslice[0]),3))
    # pyloricavgs = np.ma.masked_where(maskingslice<0.3,avgs)
    pyloricavgs = avgs[np.where(pyloricslice>=.3)][:,0::2]

    # num_nonpyloric = len(np.where(maskingslice < 0.3)[0])
    # nonpyloricavgs = np.ma.masked_where(maskingslice>=0.3,avgs)
    nonpyloricavgs = avgs[np.where(pyloricslice<.3)][:,0::2]

    if plot and (ax==None):
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.scatter(nonpyloricavgs[:,0],nonpyloricavgs[:,1],color='r',alpha=.2,s=1)
        plt.scatter(pyloricavgs[:,0],pyloricavgs[:,1],color='g',alpha=.2,s=1)
        plt.show()
    if plot and (ax!=None):
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.scatter(nonpyloricavgs[:,0],nonpyloricavgs[:,1],color='r',alpha=.2,s=1)
        ax.scatter(pyloricavgs[:,0],pyloricavgs[:,1],color='g',alpha=.2,s=1)

    solvable = False

    if evaluate:
        if high_res == False: #then, conceivable to run through and check every point in the pyloric list
            for i in range(len(pyloricavgs)):
                if progress:
                    print(i)
                dists = np.linalg.norm(nonpyloricavgs-pyloricavgs[i],axis=1)
                if np.min(dists)>tolerance:
                    solvable = True
                    if progress:
                        print(pyloricavgs[i],nonpyloricavgs[np.where(dists==np.min(dists))],np.min(dists))
                    break
        else: #takes too long to check every point, so we have a more complicated scheme (see markdown)
            elimination_pyl_list = np.copy(pyloricavgs)
            while len(elimination_pyl_list>0):
                pyl_pt = elimination_pyl_list[0]
                dists_pyl = np.linalg.norm(nonpyloricavgs-pyl_pt,axis=1)
                if np.min(dists_pyl)>tolerance:
                    solvable = True
                    if progress:
                        print(pyl_pt,nonpyloricavgs[np.where(dists_pyl==np.min(dists_pyl))],np.min(dists))
                    break
                else:
                    nonpyl_pt = nonpyloricavgs[np.where(dists_pyl==np.min(dists_pyl))]
                    dists_nonpyl = np.linalg.norm(elimination_pyl_list-nonpyl_pt,axis=1)
                    elimination_pyl_list = elimination_pyl_list[dists_nonpyl>tolerance]
                    if progress:
                        print(len(elimination_pyl_list))
            # print(solvable)
    return solvable
# print("functions defined")

recovery = np.zeros((100,5,121,5))
final_recovery = np.zeros_like(recovery[:,:,:,0])
evolved_pyl_recovery = np.zeros((100,5))
zerorange_pyl_recovery = np.zeros((100))
for i in range(100):
    for j in range(5):
        # print(np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/%s/recoverytest.dat"%(i,j))[:,2:])
        recovery[i,j] = np.loadtxt("./Specifically Evolved HP mechanisms/Every Circuit/%s/%s/recoverytest.dat"%(i,j))[:,2:]
        final_recovery[i,j] = recovery[i,j,:,-1]
        evolved_pyl_recovery[i,j] = np.sum(final_recovery[i,j]>=.3)
    hpparslice = np.loadtxt('./Specifically Evolved HP mechanisms/Every Circuit/%s/HPparslice_newrho_res5.dat'%i).T
    zerorange_pyl_recovery[i] = np.max(hpparslice)

# print("recovery data gathered")

for i in range(78,100):
    solvability = uniquevals(i,evaluate=True,plot=False)
    with open('./solvability_output.txt', 'a') as f:
        print(i,file=f)
        print("Solvable: ", solvability,file=f)
        print("Evolved solutions: ", int(np.max(evolved_pyl_recovery[i])),file=f)
        print("Zero Range Solutions: ", int(zerorange_pyl_recovery[i]),file=f)
    with open('./solvability_readable.txt', 'a') as f:
        print(i, int(solvability), int(np.max(evolved_pyl_recovery[i])), int(zerorange_pyl_recovery[i]),file=f) 