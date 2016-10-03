from particlecontainer import ParticleContainer, testfun_verify_verlet_list
from fcc import FCC
import numpy as np
import os
import pickle
import math
from logger import Log

verify = True
verify = False


def setup(nb):    
    # nb = number of unit_cells in each direction
    n_atoms = 4*nb**3 # fcc has 4 atoms per unit cell
    OK = False
    atoms_filename = 'atoms(nb={}).pickled'.format(nb)
    
    if os.path.exists(atoms_filename):
        #read the atoms from a file
        f = open(atoms_filename,'rb')
        atoms = pickle.load(f)

        if atoms.size()==n_atoms:
            OK = True 
            Log.Print('read from file:',atoms_filename)
            Log.Print("    nb =",atoms.extra['fcc_nb'])
            Log.Print("natoms =",atoms.size())
            Log.Print("  rmin =",atoms.extra['rmin'])
            Log.Print("     a =",atoms.extra['fcc_a'])
            Log.Print("cutoff =",atoms.extra['rcutoff'])
            Log.Print("     w =",atoms.extra['rcutoff'])
            Log.Print("nb*a/w =",atoms.extra['fcc_nb']*atoms.extra['fcc_a']/atoms.extra['rcutoff'])

    if not OK:
        #construct atoms and save to file
        # lattice constant
        rmin = math.pow(2,1/6)
        a = math.sqrt(2)*rmin
        offset = np.ones((3,))*(0.25*a)
        #generate atoms
        
        atoms = ParticleContainer()
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,a=a,offset=offset,verbose=False)
    
        rcutoff = 3*rmin
        
        Log.Print("    nb =",nb)
        Log.Print("natoms =",n_atoms)
        Log.Print("  rmin =",rmin)
        Log.Print("     a =",a)
        Log.Print("cutoff =",rcutoff)
        Log.Print("     w =",rcutoff)
        Log.Print("nb*a/w =",nb*a/rcutoff)

        atoms.extra = {'fcc_nb' :nb
                      ,'fcc_a'  :a
                      ,'rmin'   :rmin
                      ,'rcutoff':rcutoff
                      }
        f = open(atoms_filename,'wb')
        pickle.dump(atoms,f)
        Log.Print('written to file:',atoms_filename)
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
    if unit=='KB':
        return bytes/1024        
    if unit=='MB':
        return bytes/(1024*1024)
    if unit=='GB':
        return bytes/(1024*1024*1024)
    raise NotImplemented('unit {} not known'.format(unit))

def randomize(atoms):
    n = atoms.size()
    I = np.random.permutation(n)
    atoms.rx = atoms.rx[I]
    atoms.ry = atoms.ry[I]
    atoms.rz = atoms.rz[I]

def experiment(case,atoms,file,reserve=200,verify=False):
    """
    case : 0 = no spatial_sort
    case : 1 = no spatial_sort and atoms randomized
    case : 2 = spatial_sort
    """
    with Log("build_verlet_list -  case "+str(case)):
        rcutoff = atoms.extra['rcutoff']
        if case==0:
            if verify:
                testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve)
            else:
                atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve,verbose=False)
    
        elif case==1:
            randomize(atoms)
            if verify:
                testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve)
            else:
                atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=False,reserve=reserve,verbose=False)
            pass
            
        elif case==2:
            if verify:
                testfun_verify_verlet_list(atoms      ,rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=True ,reserve=reserve)
            else:
                atoms.build_verlet_list(algo='hilbert',rcutoff=rcutoff,cell_width=rcutoff,spatial_sort=True ,reserve=reserve,verbose=False)
    
    n_iterations = 1
    cputime = atoms.compute_interactions_verlet(verbose=True)
    while cputime<0.1:
        n_iterations = 1.5*int(math.ceil(0.1/(cputime/n_iterations)))
        cputime = atoms.compute_interactions_verlet(verbose=True,n_iterations=n_iterations)
    
    kb = memory_used(atoms,unit='KB')
    Log.Print('memory used   :',kb,'KB')
    Log.Print('total cputime :',cputime,'s')
    Log.Print('cputime/iter  :',cputime/n_iterations,'s')
    n_pairs = atoms.verlet_list.n_pairs()
    interactions_per_second = n_iterations*n_pairs/cputime
    Log.Print('interactions/s:',interactions_per_second)
    if not file is None:
        file.write(str(nb)                     +' '
                  +str(n_iterations)           +' '
                  +str(n_pairs)                +' '
                  +str(kb)                     +' '
                  +str(interactions_per_second)+'\n'
                  )

def GB(nb):
    """
    GB of memory used as a function of the number of unit cells (nb) in x/y/z direction 
    """
    return (4*6*8*nb**3)/(1024)

def one_experiment(nb,cases,verify=False):
    mem = GB(nb)
    if isinstance(cases,int):
        cases = [cases] 
    for case in cases:
        with Log('one_experiment: nb={}, mem={:.3g} kB, case={}'.format(nb,mem,case)):
            atoms = setup(nb)
            experiment(case,atoms,None,verify=verify)


if __name__=="__main__":
    one_experiment(6,range(3),verify=True)
    exit()
    
    filenames = ['data_unsorted.txt','data_randomized.txt','data_spatial_sort.txt']
    for iexp in range(3):
        with Log("experiment.py - case "+str(iexp)): 
            f = open(filenames[iexp],'a')
            nb0=3
            nb =3
            mem = GB(nb)
            while mem<0.00005*1024*1024:
                with Log('nb={}, mem={:.3g} kB'.format(nb,mem)):
                    atoms = setup(nb)
                    experiment(iexp,atoms,f)
                    
                nb = int(nb*1.26)
                if nb==nb0:
                    nb+=1
                    nb0=nb
                mem = GB(nb)
                Log.Print('')
            f.close()