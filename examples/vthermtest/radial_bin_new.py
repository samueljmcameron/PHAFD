import numpy as np


def radial_bin(nvals,N,L,dataname,nbins=20,timesteps=[0],separate=True):
    """

    Binning a dataset in fourier space.

    Params:
    -------
    
    N : int
        number of grid points in each direction assuming full fourier space
    L : double
        size of periodic domain in real space
    nproc : int
        number of processors the data is split over
    dataname : string
        name of files that contain data. Put '%' in place of processor number, and
        '*' in place of timestep. So if files are over eight processors and 30 timesteps,
        and e.g. data from processor one from timestep 4 is called dum_p1_4.txt, then
        dataname should be entered as dum_p%_*.txt
    nbins: int (optional - default is 20)
        number of bins to split data into
    timesteps : array (optional - default is [0])
        names of timesteps in your data name. If your data is only on one timestep and
        doesn't include the name, simply put in dum_p%.txt and leave this argument to
        its default value.
    separate : bool (optional - default is True)
        true if the data is split over all the processors, false otherwise
    

    Returns:
    --------

    qbins : 1D numpy array of doubles
        bin coordinates in fourier space
    flucts : 2D numpy array of doubles - len(timesteps) # of rows, len(nbins) # of cols
        binned data
    counts : 2D numpy array of ints - len(timesteps) # of rows, len(nbins) # of cols
        counts of data per bin
    

    """
    
    nx = N//2+1
    ny = N
    nz = N



    dq = 2*np.pi/L
    
    nbins = 20
    qbins = np.linspace(0,dq*N/2,num=nbins,endpoint=False)
    dqbin = qbins[1]-qbins[0]


    # data loaded in goes along x axis, then y axis, then z axis

    flucts = np.zeros([nvals,len(timesteps),nbins],float)
    counts = np.zeros([nvals,len(timesteps),nbins],int)
    
    for ti,timestep in enumerate(timesteps):


        fname = dataname.replace("*",str(timestep))
        print(fname)
        fulldata = np.loadtxt(fname,skiprows=1,delimiter=',')


        for nval in range(nvals):
            inp = fulldata[:,nval].reshape(nz,ny,nx)
            print(inp.shape)

            output = np.zeros(inp.shape)

            # shift y and z data appropriately, as the current range is
            # qx \in [0,N//2 + 1], qy \in [0,N), and qz \in [0,N), but
            # need qy \in (-N//2,N//2] and qz \in (-N//2,N//2]


            output[:N//2-1,:N//2-1,:] = inp[N//2+1:,N//2+1:,:]
            output[:N//2-1,N//2-1:,:] = inp[N//2+1:,:N//2+1,:]
            output[N//2-1:,:N//2-1,:] = inp[:N//2+1,N//2+1:,:]
            output[N//2-1:,N//2-1:,:] = inp[:N//2+1,:N//2+1,:]
        
            for iz in range(N):
                for iy in range(N):
                    for ix in range(nx):

                        q = np.sqrt(dq**2*(ix**2+(iy-N//2+1)**2+(iz-N//2+1)**2))


                        ibin = int(q/dqbin)

                        if ibin < nbins:
                            flucts[nval,ti,ibin] += output[iz,iy,ix]
                            counts[nval,ti,ibin] += 1


    return qbins+0.5*dqbin, flucts,counts
