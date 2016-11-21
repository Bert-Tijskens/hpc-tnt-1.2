import numpy as np
import matplotlib.pyplot as plt
import math


data  = np.loadtxt('hopper/experiment-simd.py.txt')
nrows = data.shape[0]//3
case0 = data[      0:  nrows,:]
case1 = data[  nrows:2*nrows,:]
case2 = data[2*nrows:3*nrows,:]

def mem(case, add_vl=False):
    """
    memory used by non contiguously accessed arrays (positions and accelerations)
    """
    n = case[:,3]
    if add_vl:
        return case[:,5]
    else:
        return n*6*8/1024

L1  = math.log10( 32)/np.log10(2)
L2 = math.log10(256)/np.log10(2)
L3= math.log10( 25*1024)/np.log10(2)

n=90000000
print(mem(case0))
ymax = 7
plt.plot( [L1,L1],[0,ymax],'-',[L2,L2],[0,ymax],'-',[L3,L3],[0,ymax],'-'
        , np.log10(mem(case0))/np.log10(2),n/case0[:,6],'o-'
        )
plt.legend(['L1','L2','L3'
           ,'case 1 (FCC regular)'],loc='upper left')
plt.xlabel('log2(memory used [kB])')
plt.ylabel('cputime relative to baseline')

plt.figure()
plt.plot( [L1,L1],[0,ymax],'-',[L2,L2],[0,ymax],'-',[L3,L3],[0,ymax],'-'
        , np.log10(mem(case0))/np.log10(2),n/case0[:,6],'o-'
        , np.log10(mem(case1))/np.log10(2),n/case1[:,6],'o-'

        )
plt.legend(['L1','L2','L3'
           ,'case 1 (FCC regular)'
           ,'case 2 (FCC permuted)'],loc='upper left')
plt.xlabel('log2(memory used [kB])')
plt.ylabel('cputime relative to baseline')

plt.figure()
plt.plot( [L1,L1],[0,ymax],'-',[L2,L2],[0,ymax],'-',[L3,L3],[0,ymax],'-'
        , np.log10(mem(case0))/np.log10(2),n/case0[:,6],'o-'
        , np.log10(mem(case1))/np.log10(2),n/case1[:,6],'o-'
        , np.log10(mem(case2))/np.log10(2),n/case2[:,6],'o-'
        )

plt.legend(['L1','L2','L3'
           ,'case 1 (FCC regular)'
           ,'case 2 (FCC permuted)'
           ,'case 3 (spatial sort)
           '],loc='upper left')
plt.xlabel('log2(memory used [kB])')
plt.ylabel('cputime relative to baseline')
plt.show()