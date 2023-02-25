import numpy as np
from genrand import polymer_coords
import matplotlib.pyplot as plt

proc = 0


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

pmer_atoms = np.array([20,30],int)

total_atoms = np.sum(pmer_atoms)

with open(f"simple_p{proc}.data","w") as writer:


    writer.write("particle_style polymer\n\n")

    writer.write(f"total_atoms {total_atoms}\n\n")


    bondlength = 0.6
    kappa = 40.0
    temp = 4.114

    Lmax = 0
    
    running_total_atoms = 0
    
    for molID,natoms in enumerate(pmer_atoms):

        L = bondlength*(natoms-1)
        if (L > Lmax):
            Lmax = L
    
        xs,ys,zs = polymer_coords(bondlength,kappa/temp,natoms,0,0,0)

        types = np.zeros([natoms],int)+1
        types[0] = types[1] = 0
        types[-2] = types[-1] = 0

        for i in range(natoms):
            writer.write(f"{i+running_total_atoms} {molID} {types[i]} {xs[i]} {ys[i]} {zs[i]} \n")

        running_total_atoms += natoms
        writer.write("\n\n")

        ax.scatter(xs,ys,zs)


ax.set_xlim(-Lmax/2,Lmax/2)
ax.set_ylim(-Lmax/2,Lmax/2)
ax.set_zlim(-Lmax/2,Lmax/2)

plt.show()
