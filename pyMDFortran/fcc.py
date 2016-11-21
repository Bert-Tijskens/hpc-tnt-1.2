import math
import numpy as np
from logger import Log

# FCC lattice points
A = np.array([[0.0, 0.0, 0.0]
             ,[0.5, 0.5, 0.0]
             ,[0.5, 0.0, 0.5]
             ,[0.0, 0.5, 0.5]])

def FCC(n,offset=[0,0,0],a=1,verbose=False):
    """
    Generate n points on a FCC lattice with a cubic lattice vector a.
    """ 
    # coordinate arrays
    x = np.zeros((n,))
    y = np.zeros((n,))
    z = np.zeros((n,))
    # translation vector
    b = np.zeros((3,))
    # xyz coordinates of i-th point
    r = np.zeros((3,)) 
    # number of unit cells that is repeated in each direction
    nb = int(math.ceil((float(n)/4)**(1./3)))
    i = 0
    for iz in range(nb):
        if i==n:
            break
        b[2] = iz + offset[2]
        for iy in range(nb):
            if i==n:
                break
            b[1] = iy + offset[1]
            for ix in range(nb):
                if i==n:
                    break
                b[0] = ix + offset[0]
                for j in range(4):
                    if i==n:
                        break
                    r = a*(b+A[j])
                    x[i] = r[0]
                    y[i] = r[1]
                    z[i] = r[2]
                    
                    i+=1
   
    if verbose:
        Log.Output('FCC:')
        Log.Output('  number of points generated:'+str(n))
        Log.Output('  lattice constant a        :'+str(a))
        Log.Output('  offset                    :'+str(offset))
        Log.Output('  maximum number of cells in each direction:'+str(nb))
        Log.Output('  x :'+str(x))
        Log.Output('  y :'+str(y))
        Log.Output('  z :'+str(z))
        Log.Output('(call fcc with verbose=False to suppress this output)')
    return x,y,z,nb

#--------------------------------------------------------------------------------------------------
# test code below
#--------------------------------------------------------------------------------------------------
if __name__=="__main__":
    x,y,z,nb = FCC(4,verbose=True)
    print(A[:,0])
    assert np.all(x==A[:,0])
    assert np.all(y==A[:,1])
    assert np.all(z==A[:,2])
    assert nb==1

    x,y,z,nb = FCC(5,verbose=True)
    assert np.all(x[:4]==A[:4,0])
    assert np.all(y[:4]==A[:4,1])
    assert np.all(z[:4]==A[:4,2])
    assert np.all(x[4]==1)
    assert np.all(y[4]==0)
    assert np.all(z[4]==0)
    assert nb==2
    
    x,y,z,nb = FCC(4,offset=[0.25,0.25,0.25],verbose=True)
    B = A + np.array([0.25,0.25,0.25])
    print(B)
    assert np.all(x==B[:,0])
    assert np.all(y==B[:,1])
    assert np.all(z==B[:,2])
    assert nb==1
    
    print('done')