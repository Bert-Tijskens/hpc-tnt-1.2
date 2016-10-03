import numpy as np
import pyHilbertCpp as hilbert 
from ParticleContainer import ParticleContainer, testfun_verify_verlet_list
from fcc import FCC
import math

"""
3374 255 1 2 4
3375 255 1 2 4
3376 255 2 2 3
3377 255 2 2 3
3378 255 2 2 4
3379 255 2 2 4
3380 255 2 2 3
3381 255 2 2 3
3382 255 2 2 4
"""
print(hilbert.ijk2h_1(1,2,4))
print(hilbert.ijk2h_1(2,2,3))
exit()
nB = 10
n_atoms = 4*nB**3
atoms = ParticleContainer(n_atoms,arrays='raih')

rmin = math.pow(2,1/6)
a = math.sqrt(2)*rmin
offset = np.ones((3,))*(0.25*a)
atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,a=a,offset=offset,verbose=False)
assert nb==nB

rcutoff = 3*rmin
atoms.compute_ijk_h(cell_width=rcutoff)
atoms.build_tables(rcutoff=rcutoff,reserve_atoms_per_cell=80,reserve_verlet=200,verbose=True)

for i in range(n_atoms):
    print(i,atoms.h[i],atoms.i[i],atoms.j[i],atoms.k[i])

H = atoms.hilbert_list.shape
for i in range(H[1]):
    print(i, atoms.hilbert_list[0,i], atoms.hilbert_list[1,i])