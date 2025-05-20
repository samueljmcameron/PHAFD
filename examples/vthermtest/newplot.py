import matplotlib.pyplot as plt
import pickle
import numpy as np


figsize =  plt.rcParams['figure.figsize']
figurename = "fluctuations"
eof=".pdf"

def noise(q,visc,dt,L,kbT=4.114):

    return 2*kbT*q**2*L**3*visc*dt



L = 256.0
visc = 0.89e-3
N = 64
dt = 1e-4
#markers = ['o','s','d','^']

with open('hists.pkl','rb') as fp:
    data = pickle.load(fp)






qbins = data[0]

for ni,(key,val) in enumerate(data[1].items()):
    print(key,val.shape)
    fig,ax = plt.subplots(figsize=figsize)
    
    ax.plot(qbins[1:],L**6*val[1:],'o',label=f"{key}")

    if key in ['c_Zx','c_Zy','c_Zz']:
        ax.plot(qbins[1:],noise(qbins[1:],visc,dt,L),'k-')

    plt.show()
