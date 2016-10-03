import numpy as np
import imp
import pyHilbert as hilbert
imp.reload(hilbert)

if __name__=="__main__":
    x = np.array([0,0,0,0,1,1,1,1],dtype=float)
    y = np.array([0,0,1,1,0,0,1,1],dtype=float)
    z = np.array([0,1,0,1,0,1,0,1],dtype=float)
    w = 1.0
    h = np.zeros(x.shape,dtype=np.int64)
    x += 0.5
    y += 0.5
    z += 0.5
    hilbert.xyzw2h(x,y,z,w,h)
    print("h:",h)
    I = np.zeros(x.shape,dtype=np.uint32)
    hilbert.sort(h,I)
    print("h:",h)
    print("I:",I)
    told = np.array([0,7,1,6,3,4,2,5],dtype=float)
    tnew = np.zeros(x.shape          ,dtype=float)
    print('told',told)
    hilbert.reorder(I,told,tnew)
    print('tnew',tnew)
    print("finished.")