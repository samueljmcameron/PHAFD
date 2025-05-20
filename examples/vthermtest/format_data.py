import sys
from radial_bin_new import radial_bin

import numpy as np
import pickle

Ns = [64]
L = 256.0
output = {}

observables = ["c_vx","c_vy","c_vz","c_Zx","c_Zy","c_Zz"]

for ni,N in enumerate(Ns):

    nprocs = 4
    timesteps = np.array(list(range(0,2)))

    dname = "converted/converted_*.csv"


    # this correction is necessary because the version of my code had an error in
    # the prefactor scaling!
    #correction = (N//2+1)*N*(N//nprocs)/N**6*L**6
        
    qbins,flucts,counts = radial_bin(6,N,L,dname,timesteps=timesteps,separate=False)


    totalflucts = {}
    totalcounts = {}
    
    for nval, name in enumerate(observables):
        totalflucts[name] = np.zeros([flucts.shape[-1]])
        totalcounts[name] = np.zeros([counts.shape[-1]],dtype=int)
        for timestep in timesteps[1:]:
            index = timestep-timesteps[0]

            totalflucts[name] += flucts[nval,index,:]
            totalcounts[name] += counts[nval,index,:]

        output[name] = np.where(totalcounts[name] > 0,
                                totalflucts[name]/totalcounts[name],0)


with open('hists.pkl', 'wb') as fp:
    pickle.dump([qbins,output],fp)
