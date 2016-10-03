import numpy as np
import math

N = 100

def lj2(r2):
    r2i=1./r2
    r6i=r2i*r2i*r2i
    lj = 4.*r6i*(r6i-1.)
    return lj

rm = math.pow(2,1./6)
rm2 = rm**2

r2 = np.zeros((N,),dtype=np.float32)
for i in range(N):
    r2[i]=i*i*rm2

for i in range(1,N):
    LJ = np.float32(lj2(r2[i]))
    print("{} {:20.15f}".format(i,LJ))

