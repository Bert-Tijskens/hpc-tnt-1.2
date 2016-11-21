from __future__ import division
import particlecontainer as pc
from fcc import FCC
import numpy as np
import os
import pickle
import math
from logger import Log

verify = True
verify = False

def setup(nb,case,reserve=200,verify=False,force_construct=False):    
    # nb = number of unit_cells in each direction
    n_atoms = 4*nb**3 # fcc has 4 atoms per unit cell
    OK = False
    atoms_filename = 'atoms(nb={},case={}).pickled'.format(nb,case)
    
    if os.path.exists(atoms_filename) and not force_construct:
        #read the atoms from a file
        f = open(atoms_filename,'rb')
        atoms = pickle.load(f)

        if atoms.size()==n_atoms:
            OK = True 
            Log.Output('read from file:',atoms_filename)
            Log.Output("    nb =",atoms.extra['fcc_nb'])
            Log.Output("natoms =",atoms.size())
            Log.Output("  rmin =",atoms.extra['rmin'])
            Log.Output("     a =",atoms.extra['fcc_a'])
            Log.Output("cutoff =",atoms.extra['rcutoff'])
            Log.Output("     w =",atoms.extra['rcutoff'])
            Log.Output("nb*a/w =",atoms.extra['fcc_nb']*atoms.extra['fcc_a']/atoms.extra['rcutoff'])

    if not OK:
        #construct atoms and save to file
        # lattice constant
        rmin = math.pow(2,1/6)
        a = math.sqrt(2)*rmin
        offset = np.ones((3,))*(0.25*a)
        #generate atoms
        
        atoms = pc.ParticleContainer()
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,a=a,offset=offset,verbose=False)
    
        rcutoff = 3*rmin
        
        Log.Output("    nb =",nb)
        Log.Output("natoms =",n_atoms)
        Log.Output("  rmin =",rmin)
        Log.Output("     a =",a)
        Log.Output("cutoff =",rcutoff)
        Log.Output("     w =",rcutoff)
        Log.Output("nb*a/w =",nb*a/rcutoff)

        atoms.extra = {'fcc_nb' :nb
                      ,'fcc_a'  :a
                      ,'rmin'   :rmin
                      ,'rcutoff':rcutoff
                      }
        
        with Log("build_verlet_list(atoms={}, case={})".format(n_atoms,case)):
            rcutoff = atoms.extra['rcutoff']
            if case==0:
                if verify:
                    pc.testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve)
                else:
                    atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve,verbose=False)
        
            elif case==1:
                randomize(atoms)
                if verify:
                    pc.testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve)
                else:
                    atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve,verbose=False)
                pass
                
            elif case==2:
                if verify:
                    pc.testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=True ,reserve=reserve)
                else:
                    atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=True ,reserve=reserve,verbose=False)

        atoms.extra['mem_kB'] = memory_used(atoms,unit='kB')
                    
#         f = open(atoms_filename,'wb')
#         pickle.dump(atoms,f)
#         Log.Output('written to file:',atoms_filename)
    return atoms
    
def memory_used(atoms,precision=2,unit="B"):
    """
    return memory used in Bytes for the verlet list case
    """
    n_atoms = atoms.size()
    n_arrays = 6
    bytes = n_atoms*n_arrays*4*precision
    n_pairs = atoms.verlet_list.n_pairs()
    bytes += (n_pairs+n_atoms)*4
    if unit =='B':
        return bytes
    if unit=='kB':
        return bytes/1024        
    if unit=='MB':
        return bytes/(1024*1024)
    if unit=='GB':
        return bytes/(1024*1024*1024)
    raise NotImplementedError('unit {} not known'.format(unit))

def randomize(atoms):
    n = atoms.size()
    I = np.random.permutation(n)
    atoms.rx = atoms.rx[I]
    atoms.ry = atoms.ry[I]
    atoms.rz = atoms.rz[I]

def experiment(case,atoms,file):
    """
    case : 0 = no spatial_sort
    case : 1 = no spatial_sort and atoms randomized
    case : 2 = spatial_sort
    """
    
    n_iterations = 1
    cputime = atoms.compute_interactions(verbose=True)
    cputime_desired = 10. #seconds 
    while cputime<cputime_desired:
        if cputime==0:
            n_iterations *= 100
        else:
            n_iterations = 1.5*int(math.ceil(cputime_desired/(cputime/n_iterations)))
        cputime = atoms.compute_interactions(verbose=True,n_iterations=n_iterations)
    
    kb = memory_used(atoms,unit='kB')
    Log.Output('memory used   :',kb,'kB')
    Log.Output('total cputime :',cputime,'s')
    Log.Output('cputime/iter  :',cputime/n_iterations,'s')
    n_atoms = atoms.size()
    n_pairs = atoms.verlet_list.n_pairs()
    interactions_per_second = n_iterations*n_pairs/cputime
    Log.Output('interactions/s:',interactions_per_second)
    if not file is None:
        file.write(str(case)                    +' '
                  +str(atoms.extra['fcc_nb'])   +' '
                  +str(n_iterations)            +' '
                  +str(n_atoms)                 +' '
                  +str(n_pairs)                 +' '
                  +str(kb)                      +' '
                  +str(interactions_per_second) +'\n'
                  )
        file.flush()

def kB(nb):
    """
    kB of memory used as a function of the number of unit cells (nb) in x/y/z direction 
    """
    return (4*6*8*nb**3)/(1024)

def one_experiment(nb,cases,verify=False,file=None):
    if isinstance(cases,int):
        cases = [cases] 
    for case in cases:
        with Log('one_experiment: nb={}, mem={:.3g} kB, case={}'.format(nb,kB(nb),case)):
            atoms = setup(nb,case,verify=verify)
            experiment(case,atoms,file)

def vl_timing(nb):
    use_cpp = pc.use_cpp
    
    atoms = setup(nb,case=2,force_construct=True) # with spatial sort
    pc.use_cpp = False
    atoms = setup(nb,case=2,force_construct=True) # with spatial sort
    
    pc.use_cpp = use_cpp
    
if __name__=="__main__":
#     vl_timing(12)
#     exit()
    f = open('experiment-simd.py.txt','w')
    
#     one_experiment(81,2,verify=False,file=f)
#     exit()

            
    with Log("experiment.py",out='experiment.py.log'):
        for case in range(3):
            nb0=3
            nb =3
            mem = kB(nb)
            while mem<0.7*1024*1024: # * 1kB
                with Log('experiment: nb={}, mem={:.3g} kB, case={}'.format(nb,mem,case)): 
                    atoms = setup(nb,case)
                    experiment(case,atoms,f)
                    
                nb = int(nb*1.26)
                if nb==nb0:
                    nb+=1
                    nb0=nb
                mem = kB(nb)
                Log.Output('')
    f.close()