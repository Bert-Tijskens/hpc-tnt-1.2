from pyMDFortran import md as md
import random
import numpy as np
 
if __name__=='__main__':
    i0 = 9
    i1 = 29
    nk1 = 2**(i1-1)
    print(nk1, nk1*6*8/(1024**3), 'GB')
    rx = np.random.random((nk1,))
    ry = np.random.random((nk1,))
    rz = np.random.random((nk1,))
    ax = np.zeros ((nk1,),dtype=np.float64)
    ay = np.zeros ((nk1,),dtype=np.float64)
    az = np.zeros ((nk1,),dtype=np.float64)
    f = open('ppmd01ff.py.txt','w')
    for i in range(i0,i1):
        n = 2**i
        k = (2**(i1-1))//n
        nk = n*k
        assert nk==nk1
        print(i)
        f.write("{} {} {} {} {} {}\n".format(i,n,k,n*k, 6*n*8/1024,md.ppmd01ff(rx,ry,rz,ax,ay,az,n,k)))

    print('done')