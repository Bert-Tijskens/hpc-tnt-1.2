from pymd import mymd as md


import numpy as np
import time
 

xyz0 = np.random.random( (3,) )



def uniform_real(n,a=0.,b=1.,dtype=np.float):
    r = np.random.random( (n,) )
    if a!=0:
        r = r*(b-a) 
        r += a
    elif b!=1:
        r = r*b  
#     print(r)
    if dtype!=r.dtype:
        rc = r.astype(dtype)
#         print('casting')
        return rc
    else:
        return r      

def experiment(m,k):
    dtype = np.float64
    rx = uniform_real(m,1,3,dtype=dtype)
    ry = uniform_real(m,1,3,dtype=dtype)
    rz = uniform_real(m,1,3,dtype=dtype)
    start = time.time()
#     print(rx.__array_interface__['data'],type(rx.__array_interface__['data'][0]))
    for ik in range(k):
        e = md.interaction_energy(rx,ry,rz,xyz0[0],xyz0[1],xyz0[2])
#         print(ik,e)
        
    stop = time.time()
    print(m,k,"walltime:",stop-start)

if __name__=="__main__":
    n = 20
    m = 2**n
    k =1
    for i in range(1,n-8):
        experiment(m,k)
        m=m/2
        k=k*2
    
    print("done")