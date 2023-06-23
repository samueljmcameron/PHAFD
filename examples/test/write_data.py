import sys
import numpy as np
path = sys.argv[1] + '/Documents/polymerode/scripts/'
sys.path.append(path)
import subprocess
from set_params import set_params
from convertVTKtoNParrays import convertVTKtoNParrays
from fix_roundoff import fix_roundoff
from genpoints import genpoints
## first generate an equilibrium (ish) configuration from
## running the executable

P = 50.0
bondlength = 4.0
temp = 4.114

N = int(sys.argv[2])
seed = sys.argv[3]
dt = sys.argv[4]
run = sys.argv[5]
attract = sys.argv[6]
proc = int(sys.argv[7])



npolymers = 1


folder = f'polymercoords_{attract}'

L = (N-1)*bondlength
Re2e = np.sqrt(2*L*P*(1-P/L*(1-np.exp(-L/P))))

if (proc < npolymers):
    tparams = set_params(N,bondlength=bondlength,kbT=temp)



    x1 = np.array([-40.0,-40.0,-40.0])
    xN = -1*x1


    fname = folder + f'/run{run}_pmer{proc}'
    args= f"{sys.argv[1]}/poly/bin/beadrodpmer-bin -in initpolymers.in -var beads {N} -var bondlength {bondlength} -var zeta_para {tparams['friction_para']} -var zeta_perp {tparams['friction_perp']} -var x1x {x1[0]} -var x1y {x1[1]} -var x1z {x1[2]} -var xNx {xN[0]} -var xNy {xN[1]} -var xNz {xN[2]} -var bending {P*temp} -var temp {temp} -var seed {seed} -var dt {dt} -var fname {fname}".split()
    
    popen = subprocess.Popen(args,stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    #print(output.decode('UTF-8'))

    xs,ys,zs = convertVTKtoNParrays(f'{fname}_0.vtp')
    xs,ys,zs = fix_roundoff(xs,ys,zs)

    nbeads = xs.size
    
    mols = np.zeros([nbeads],int)+proc
    types = np.copy(mols)*0 + 1

    types[::20] = 0
    types[1::20] = 0
    types[2::20] = 0    
    types[3::20] = 0    
    types[4::20] = 0    
    types[5::20] = 0    
    types[6::20] = 0    
    types[7::20] = 0    
    types[8::20] = 0    
    types[9::20] = 0    

    with open(f'{folder}/run{run}_empty_p{proc}.data','w') as myfile:

        for i in range(nbeads):
            ID = i + proc*nbeads
            myfile.write(f"{ID} {mols[i]} {types[i]} {xs[i]} {ys[i]} {zs[i]}\n")

else: # just write empty file
    with open(f'{folder}/run{run}_empty_p{proc}.data','w') as myfile:
        myfile.write(f"\n")

if (proc == 0):
    inputargs= f"-var bondlength {bondlength} -var zeta_para {tparams['friction_para']} -var zeta_perp {tparams['friction_perp']} -var bending {P*temp} -var temp {temp} -var ljcut {bondlength} -var ljeps {1.0} -var ljsigma {bondlength/2**(1./6.)} -var nskin {bondlength/2.0} -var dt {dt}"

    sys.stdout.write(inputargs)
