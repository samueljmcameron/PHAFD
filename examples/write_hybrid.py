import numpy as np


nprocs = 4



L = 100

allPmerNs = [2,1,1,1]
allSphereNs = [4,12,12,12]


allAtomsPerPmer = [[12,30],[20],[18],[30]]

allPmerIDs = [["0-12","12-42"],["42-62"],["62-80"],["80-110"]]

allSphereIDs = [list(range(110,114)),list(range(114,126)),
                list(range(126,138)),list(range(138,150))]


allPmerMolIDs = [[0,1],[2],[3],[4]]

allSphereMolIDs = [list(range(5,9)),list(range(9,21)),
                   list(range(21,33)),list(range(33,45))]

allPmerFiles = [["polymercoords/first_nt_0.vtp","polymercoords/first_dt_0.vtp"],
                ["polymercoords/second_nt_0.vtp"],["polymercoords/first_st_0.vtp"],
                ["polymercoords/second_st_0.vtp"]]

allPmerTypes = [["0*2 1*8 0*2","1*14 0*2 1*14"],["1*18 0*2"],["1*18"],["1*20 0*10"]]

allSphereTypes = [[2]*4,[2]*12,[2]*12,[2]*12]




for proc in range(nprocs):


    npmers = allPmerNs[proc]
    nspheres = allSphereNs[proc]

    pmerIDs = allPmerIDs[proc]
    sphereIDs = allSphereIDs[proc]

    pmerMolIDs = allPmerMolIDs[proc]
    sphereMolIDs = allSphereMolIDs[proc]

    pmerTypes = allPmerTypes[proc]
    sphereTypes = allSphereTypes[proc]

    atomsPerPmer = allAtomsPerPmer[proc]
    pmerFiles = allPmerFiles[proc]
    
    xs = (np.random.rand(nspheres)-0.5)*L
    ys = (np.random.rand(nspheres)-0.5)*L
    zs = (np.random.rand(nspheres)-0.5)*L
    radii = np.zeros([nspheres],float)+0.5


    with open(f"hybridpython_p{proc}.data","w") as writer:

    
        writer.write("particle_style hybrid polymer sphere\n\n")

        writer.write(f"particlesperstyle {npmers} {nspheres}\n\n")

        writer.write("atomsperpolymer")

        for plen in atomsPerPmer:
            writer.write(f" {plen}")

        writer.write("\n\n")

        for pID, pMolID,pType,pFile in zip(pmerIDs,pmerMolIDs,pmerTypes,pmerFiles):
            writer.write(f"{pID} {pMolID} {pType} {pFile} \n")

        writer.write("\n\n")

        for sID,sMolID,sType,x,y,z,r in zip(sphereIDs,sphereMolIDs,sphereTypes,xs,ys,zs,radii):
            writer.write(f"{sID} {sMolID} {sType} {x} {y} {z} {r}\n")
