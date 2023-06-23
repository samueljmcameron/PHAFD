import numpy as np

def genpoints(R,a,N):
    """
    generate a discrete set of points that are equi-spaced and make a curve
    with curvate R, the centre of the curve is set along the line (1,1,1)
    and the direction of the curve is parallel to azimuthal angle = pi/4
    in standard spherical coordinates

    params:
    
    R : float
        radius of curvature for the curve
    a : float
        spacing between points in the curve
    N : integer
        number of points along the curve (should be odd)

    """

    L =(N-1)*a
    if L > 2*np.pi*R:
        raise ValueError('Length of polymer is longer than circumference of sphere')

    phi = 2*np.arcsin(0.5*a/R)

    xs = np.empty([N],float)
    ys = np.empty([N],float)
    zs = np.empty([N],float)*0


    for i in range(N):
        xs[i] = R*np.cos(i*phi)
        ys[i] = R*np.sin(i*phi)

    middle_angle = N//2*phi

    dangle = np.pi/4.-middle_angle

    
    # rotation to centre at phi = pi/4.
    R1 = np.array([[np.cos(dangle),-np.sin(dangle),0],
                   [np.sin(dangle),np.cos(dangle),0],
                   [0,0,1]])

    # rotation of basis vectors so that new x, y axes are rotated
    # by -pi/4.
    Q = np.array([[np.sqrt(2)/2.,-np.sqrt(2)/2.,0],
                  [np.sqrt(2)/2.,np.sqrt(2)/2.,0],
                  [0,0,1]])
    
    # rotation of pi/2. about new y axis
    Rtilde2 = np.array([[0,0,1.],[0,1.,0],[-1.,0,0]])

    # rotation of of np.arctan(np.sqrt(2)) about new x axis

    ang = np.pi/2.-np.arctan(np.sqrt(2))
    Rtilde3 = np.array([[1,0,0],
                        [0,np.cos(ang),-np.sin(ang)],
                        [0,np.sin(ang),np.cos(ang)]])
    
    
    # full rotation:

    rot = Q.transpose() @ Rtilde3 @ Rtilde2 @ Q @ R1
    

    
    for i in range(N):
        # rotate vectors counter clockwise by the angle so
        # that middle of curve is at phi = pi/4.
        vec = np.array([xs[i],ys[i],zs[i]])
        
        newvec = rot @ vec

        xs[i] = newvec[0]
        ys[i] = newvec[1]
        zs[i] = newvec[2]
    

    xm,ym,zm = xs[N//2],ys[N//2],zs[N//2]
    rmid = np.sqrt(xm**2+ym**2+zm**2)
    thetamid = np.arctan(np.sqrt(xm**2+ym**2)/zm)
    phimid = np.arctan(ym/xm)
    #print(rmid,thetamid,phimid)
        
    return xs,ys,zs


if __name__ == "__main__":

    bondlength = 4.0

    N = 51

    # center of sphere
    xc = np.array([7,0,0],float)
    radius=88.0

    xs,ys,zs = genpoints(radius,bondlength,N)

    xs += xc[0]
    ys += xc[1]
    zs += xc[2]

    nbeads = xs.size


    with open(f'sample.csv','w') as myfile:
        myfile.write(f"x,y,z\n")
        for i in range(nbeads):
            myfile.write(f"{xs[i]},{ys[i]},{zs[i]}\n")
