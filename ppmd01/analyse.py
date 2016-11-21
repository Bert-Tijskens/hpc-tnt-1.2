import numpy as np
import matplotlib.pyplot as plt
import math

n = 2**29
B = n*(3*8+1*4)
kB = B/1024.
MB = kB/1024.
GB = MB/1024.
print("memory used = {} GB".format(GB))

data  = np.loadtxt('ppmd01-no-xHost.txt')
datax = np.loadtxt('ppmd01-xHost.txt')

x_memory_used = True
if x_memory_used:
    L1  = math.log10( 32)/np.log10(2)
    L2 = math.log10(256)/np.log10(2)
    L3= math.log10( 25*1024)/np.log10(2)
    x = np.log10(24*(data[:,0]/1024))/np.log10(2)
else:
    L1  = math.log10( 32*1024)
    L2 = math.log10(256*1024)
    L3= math.log10(25*1024*1024)
    x = np.log10(data[:,0])

print(L1)
print(L2)
print(L3)
plt.plot( [L1,L1],[0,25],'-',[L2,L2],[0,25],'-',[L3,L3],[0,25],'-'
        , x, data[:,2], 'o-'
        , x, data[:,3], 'o-'
        , x, data[:,4], 'o-'
        , x,datax[:,4], 'o-'
        )
plt.legend(['L1','L2','L3', "contiguous SoA",'contiguous AoS','random access','random -xHost'],loc='upper left')
if x_memory_used:
    plt.xlabel('log2(mem used [kB])')
else:
    plt.xlabel('log10(array length)')
plt.ylabel('cpu time [s]')
plt.show()
exit()
