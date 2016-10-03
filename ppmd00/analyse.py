import numpy as np
import matplotlib.pyplot as plt
import math

n = 2**29
B = n*(3*8+1*4)
kB = B/1024.
MB = kB/1024.
GB = MB/1024.
print("memory used = {} GB".format(GB))

data = np.loadtxt('ppmd01a.txt')

x32  = math.log10( 32*1024)
x256 = math.log10(256*1024)
x25MB= math.log10(25*1024*1024)

plt.plot([x32,x32],[0,25],'-',[x256,x256],[0,25],'-',[x25MB,x25MB],[0,25],'-',np.log10(data[:,0]),data[:,2],'o-',np.log10(data[:,0]),data[:,3],'o-')
plt.legend(['L1','L2','L3', "contiguous access",'random access'])
plt.xlabel('log10(# LJ potential evaluations)')
plt.ylabel('cpu time [s]')
plt.show()