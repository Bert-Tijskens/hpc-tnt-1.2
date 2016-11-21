import numpy as np
import matplotlib.pyplot as plt
import math

data  = np.loadtxt('ppmd01ff.py.txt')

L1  = math.log10( 32)/np.log10(2)
L2 = math.log10(256)/np.log10(2)
L3= math.log10( 25*1024)/np.log10(2)

plt.plot( [L1,L1],[0,25],'-',[L2,L2],[0,25],'-',[L3,L3],[0,25],'-'
        , np.log10(data[:,4])/np.log10(2),data[:,5],'o-')

plt.legend(['L1','L2','L3','contiguous access'],loc='upper left')
plt.xlabel('log2(memory used [kB])')
plt.ylabel('cputime [s]')
plt.show()