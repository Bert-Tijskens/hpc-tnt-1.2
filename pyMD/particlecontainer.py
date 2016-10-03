import sys
import numpy as np
import pyHilbertCpp as hilbert 
import pyMDFortran

from logger import Log,ProgressBar
from verletlist import VerletList

import imp
imp.reload(hilbert)
imp.reload(pyMDFortran)

md = pyMDFortran.md
atom_indices_unsigned = False

E = np.array([
             #[[0,0,0],[0,0,0]],
              [[0,0,0],[1,0,0]], #0
              [[0,0,0],[0,1,0]], 
              [[0,0,0],[0,0,1]], 
              [[0,0,0],[1,1,0]], 
              [[0,0,0],[1,0,1]], 
              [[0,0,0],[0,1,1]], 
              [[0,0,0],[1,1,1]], 
              [[1,0,0],[0,1,0]], #7
              [[1,0,0],[0,0,1]], 
              [[1,0,0],[0,1,1]], 
              [[0,1,0],[0,0,1]], #10 
              [[0,1,0],[1,0,1]],
              [[0,0,1],[1,1,0]]  #12
                               ],dtype=np.int32)


class ParticleContainer():
    """
    """
    def __init__(self,n=0,real_type=np.float64,arrays='rvaih'):
        self.real_type = real_type
        self.reserve = 50
        if n==0:
            # arrays will be created on the fly... 
            # see the test code below for an example
            return
        if arrays=='all':
            arrays_ = 'rvaih'
        else:
            arrays_ = arrays
        # position coordinates
        if 'r' in arrays_:        
            self.rx = np.zeros((n,),dtype=real_type)
            self.ry = np.zeros((n,),dtype=real_type)
            self.rz = np.zeros((n,),dtype=real_type)
        # velocity coordinates
        if 'v' in arrays_:
            self.vx = np.zeros((n,),dtype=real_type)
            self.vy = np.zeros((n,),dtype=real_type)
            self.vz = np.zeros((n,),dtype=real_type)
        # acceleration coordinates
        if 'a' in arrays_:
            self.ax = np.zeros((n,),dtype=real_type)
            self.ay = np.zeros((n,),dtype=real_type)
            self.az = np.zeros((n,),dtype=real_type)
        # cell coordinates [ijk ijk ijk ...]
        if 'i' in arrays_:
            self.i = np.zeros((n,),dtype=np.int32) 
            self.j = np.zeros((n,),dtype=np.int32) 
            self.k = np.zeros((n,),dtype=np.int32) 
        # Hilbert indices
        if 'h' in arrays_:
            self.h  = np.zeros((n,),dtype=np.int64)
            
    def size(self):
        if not hasattr(self,'rx'):
            return 0
        else:
            return self.rx.size

    def resize(self,n=None,arrays=None):
        """
        """
        if n is None:
            # rx,ry,rz have been provided from outside the class
            N = self.size()
        else:
            N = n
            self.rx = np.empty((N,),dtype=self.real_type)
            self.ry = np.empty((N,),dtype=self.real_type)
            self.rz = np.empty((N,),dtype=self.real_type)
            
        if n==0:
            raise RuntimeError("resize to n=0.")
        
        if (hasattr(self,'vx') and self.vx.size!=N) or ('v' in arrays and not hasattr(self,'vx')):
            self.vx = np.empty((N,),dtype=self.real_type)
            self.vy = np.empty((N,),dtype=self.real_type)
            self.vz = np.empty((N,),dtype=self.real_type)
        if (hasattr(self,'ax') and self.ax.size!=N) or ('a' in arrays and not hasattr(self,'ax')):
            self.ax = np.empty((N,),dtype=self.real_type)
            self.ay = np.empty((N,),dtype=self.real_type)
            self.az = np.empty((N,),dtype=self.real_type)
        if (hasattr(self,'i') and self.i.size!=N) or ('i' in arrays and not hasattr(self,'i')):
            self.i = np.empty((N,),dtype=np.int32)
            self.j = np.empty((N,),dtype=np.int32)
            self.k = np.empty((N,),dtype=np.int32)
        if hasattr(self,'h') and self.h.size!=N or ('h' in arrays and not hasattr(self,'h')):
            self.h =  np.empty((N,),dtype=np.int64)
        
    def compute_ijk(self,cell_width=1):
        """
        Compute the cell indices (ijk) of all atoms for cell width cell_width>
        """
        self.cell_width = cell_width
        self.i = np.int32(np.floor(self.rx/cell_width))
        self.j = np.int32(np.floor(self.ry/cell_width))
        self.k = np.int32(np.floor(self.rz/cell_width))
        
    def compute_ijk_h(self,cell_width=1):
        """
        Compute cell indices (ijk) and Hilbert indices of all atoms for cell width cell_width.
        """
        self.cell_width = cell_width
        
        self.resize(arrays='ijkh')
        
        if self.rx.dtype==np.float32:
            hilbert.xyzw2ijkh_float32(self.rx,self.ry,self.rz,np.float32(cell_width),self.i,self.j,self.k,self.h)
        elif self.rx.dtype==np.float64:
            hilbert.xyzw2ijkh_float64(self.rx,self.ry,self.rz,np.float64(cell_width),self.i,self.j,self.k,self.h)
        else:
            raise NotImplemented('Position coordinates must have numpy.ndarray.dtype np.float32 or np.float64.')

    def compute_h(self,cell_width=1):
        """
        Compute Hilbert indices of all atoms for cell width cell_width.
        """
        self.cell_width = cell_width
        
        self.resize(arrays='h')        
        
        if self.rx.dtype==np.float32:
            hilbert.xyzw2h_float32(self.rx,self.ry,self.rz,np.float32(cell_width),self.h)
        elif self.rx.dtype==np.float64:
            hilbert.xyzw2h_float64(self.rx,self.ry,self.rz,np.float64(cell_width),self.h)
        else:
            raise NotImplemented('Position coordinates must have numpy.ndarray.dtype np.float32 or np.float64.')
        

    def spatial_sort(self):         
        """
        spatial sort of all arrays by hilbert index
        """   
        # 1. indirect sort of the hilbert indices
        self.compute_h(self.cell_width)
        I = np.argsort(self.h)
        # sort all arrays:
        self.h = self.h[I]
        if hasattr(self,'rx'):
            self.rx = self.rx[I]
            self.ry = self.ry[I]
            self.rz = self.rz[I]
        if hasattr(self,'vx'):
            self.vx = self.vx[I]
            self.vy = self.vy[I]
            self.vz = self.vz[I]
        if hasattr(self,'ax'):
            self.ax = self.ax[I]
            self.ay = self.ay[I]
            self.az = self.az[I]
        if hasattr(self,'i'):
            self.i = self.i[I]
            self.j = self.j[I]
            self.k = self.k[I]
            
    def verlet_list_add_cell_(self,h0,rcutoff2,I):
        """
        update verlet lists for h0-h0 interactions (intra-cell)
        """
        first_atom_in_h0 = self.hilbert_list[h0,0]
        n_atoms_in_h0    = self.hilbert_list[h0,1]
        for ia in range(first_atom_in_h0+1,first_atom_in_h0+n_atoms_in_h0):
            if I is None:
                i_atom = ia
            else:
                i_atom = I[ia]
            xi = self.rx[i_atom]
            yi = self.ry[i_atom]
            zi = self.rz[i_atom]
            for ja in range(first_atom_in_h0,ia):
                if I is None:
                    j_atom = ja
                else:
                    j_atom = I[ja]
                r2 = (xi-self.rx[j_atom])**2 + (yi-self.ry[j_atom])**2 + (zi-self.rz[j_atom])**2
                if r2<rcutoff2: 
                    self.verlet_list.add(i_atom,j_atom)
#                                 logger.print('{}:{}-{}'.format(ijk0,i_atom,j_atom))
        
    def verlet_list_add_cell_cell_(self,h0,h1,rcutoff2,I):
        """
        update verlet lists for h0-h1 interactions (h0!=h1)
        """
        assert h0!=h1
        first_atom_in_h0 = self.hilbert_list[h0,0]
        n_atoms_in_h0    = self.hilbert_list[h0,1]
        first_atom_in_h1 = self.hilbert_list[h1,0]
        n_atoms_in_h1    = self.hilbert_list[h1,1]
        for ia in range(first_atom_in_h0,first_atom_in_h0+n_atoms_in_h0):
            if I is None:
                i_atom = ia
            else:
                i_atom = I[ia]
            xi = self.rx[i_atom]
            yi = self.ry[i_atom]
            zi = self.rz[i_atom]
            for ja in range(first_atom_in_h1,first_atom_in_h1+n_atoms_in_h1):
                if I is None:
                    j_atom = ja
                else:
                    j_atom = I[ja]
                r2 = (xi-self.rx[j_atom])**2 + (yi-self.ry[j_atom])**2 + (zi-self.rz[j_atom])**2
                if r2<rcutoff2: 
                    self.verlet_list.add(i_atom,j_atom)
#                     added = '+'
#                 else:
#                     added = ''
#                 Log.current_logger.print('{}-{}:{}-{} : {} <? {} {}'.format(h0,h1,i_atom,j_atom,r2,rcutoff2,added))
#         Log.current_logger.print('\n')

    def set_attribute_(self,name,value):
        if value is None:
            if hasattr(self,name):
                pass
            else:
                raise RuntimeError('ParticleContainer.build_verlet_list: {} not provided.'.format(name))
        else:
            setattr(self, name, value)

    def build_verlet_list(self,algo='hilbert',rcutoff=None,reserve=50,verbose=False,cell_width=None,spatial_sort=True):
        if rcutoff is None:
            raise RuntimeError('ParticleContainer.build_verlet_list: rcutoff parameter must provided.')
            
        if algo=='hilbert':
            if cell_width is None:
                cw = rcutoff
            else:
                cw = cell_width

            if cw==rcutoff:
                cw *= 1.000000000001
            assert cw > rcutoff, 'ParticleContainer.build_verlet_list(algo="hilbert",...): The cell_width parameter must be larger than the rcutoff parameter'
                
            self.cell_width= cw
            
            return self.build_verlet_list_hilbert(self.cell_width,rcutoff,reserve,verbose,spatial_sort=spatial_sort)
        
        elif algo=='brute_force':
            return self.build_verlet_list_brute_force(rcutoff,reserve,verbose)
        
        else:
            raise NotImplementedError('algo='+str(algo))
        
    def build_verlet_list_hilbert(self,cell_width,rcutoff,reserve,verbose,spatial_sort=True):
        """
        1. perform a hilbert sort
        2. construct hilbert list
        3. construct verlet list
        """        
        with Log('ParticleContainer.build_verlet_list_hilbert: w={}, rc={}'.format(cell_width,rcutoff)):
                        
            # 1. indirect sort of the hilbert indices
            with Log("ParticleContainer.build_verlet_list - spatial sort"):
                if spatial_sort:
                    self.spatial_sort()
                    h = self.h
                    I = None
                else:
                    Log.Print('Spatial sort skipped on request.')
                    self.compute_h(self.cell_width)
                    I = np.argsort(self.h)
                    h = self.h[I]
#                     Log.Print('h',h)
#                     Log.Print('I',I)

            # now h is sorted array of hilbert indices
            # and I is the indirection array if the particles are not spatially sorted, otherwise I = None
            n_atoms = self.size()
            hend = np.max(h)+1
            # 2. construct hilbert list
            with Log("ParticleContainer.build_verlet_list - building hilbert list") as logger:
                self.hilbert_list = np.empty((hend,2),dtype=np.int32)
                #   hilbert_list[0,h] is the index of the first atom in the cell with hilbert index h, 
                #   hilbert_list[1,h] is the number of atoms in that cell
                ia = 0
                for ih in range(hend):
                    self.hilbert_list[ih,0] = ia
                    n_atoms_in_h = 0
                    while ia<n_atoms and h[ia]==ih:
                        n_atoms_in_h += 1
                        ia += 1
                    self.hilbert_list[ih,1] = n_atoms_in_h
#                 Log.Print('hilbert_list',self.hilbert_list)
                                    
            # 3. construct verlet list
            if verbose:
                progressBar = ProgressBar()
            else:
                progressBar = None
            with Log("ParticleContainer.build_verlet_list - building verlet list", progressBar=progressBar) as logger:
                rcutoff2 = rcutoff**2
                self.verlet_list = VerletList( self.size(), algo='hilbert(cw={}, rc={})'.format(cell_width,rcutoff), reserve=reserve,rcutoff=rcutoff)
                #loop over all cells
                ijk00 = np.empty((3,),dtype=np.int32)
                for h00 in range(hend):
                    hilbert.h2ijk_1(h00,ijk00)
                    h0 = h00
                    ijk0 = ijk00.copy()
                    # intra-cell
                    self.verlet_list_add_cell_(h0, rcutoff2,I)
                            
                    # loop over neighbouring cells
                    nb0 = [0,7,10,12,13]
                    for inb0 in range(4):
                        ijk0 = ijk00 + E[nb0[inb0],0,:]
                        h0 = hilbert.ijk2h_1(ijk0)
#                         logger.print("h0={}:{}".format(h0,ijk0))
                        if -1<h0<hend:
                            for nb in range(nb0[inb0],nb0[inb0+1]):
#                                 print(nb)
                                ijk1 = ijk00 + E[nb,1,:]
                                h1 = hilbert.ijk2h_1(ijk1)
    #                             logger.print("h1={}:{}".format(h1,ijk1))
                                if -1<h1<hend:
    #                                 logger.print('{}:{}-{}:{}'.format(h0,ijk0,h1,ijk1))
                                    self.verlet_list_add_cell_cell_(h0,h1,rcutoff2,I)
                    
        return self.verlet_list

    def build_verlet_list_brute_force(self, rcutoff, reserve, verbose ):
        """
        brute force O(N*N) verlet list construction.
        mainly for testing the Verlet list constructed by build_tables. 
        """
        if not verbose: p=None
        else:           p=ProgressBar()
            
        with Log('ParticleContainer.build_verlet_list_brute_force',progressBar=p) as logger:
            n_atoms = self.size()
            self.verlet_list = VerletList(n_atoms,algo='BruteForce(rc={})'.format(rcutoff),reserve=reserve,rcutoff=rcutoff)
            rcutoff2 = rcutoff*rcutoff
            for i_atom in range(n_atoms):
                logger.show_progress(i_atom/n_atoms)
                xi = self.rx[i_atom]
                yi = self.ry[i_atom]
                zi = self.rz[i_atom]
                for j_atom in range(i_atom):
                    xj = self.rx[j_atom]
                    yj = self.ry[j_atom]
                    zj = self.rz[j_atom]
                    r2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                    if r2<rcutoff2:
                        self.verlet_list.add(i_atom,j_atom)
        return self.verlet_list
                    
    def print_verlet_list_in_fortran_(self):
        """
        Print the verlet list from within the fortran module 
        To verify that the fortran module does the correct thing 
        """
        with Log('print_verlet_list_in_fortran_') as logger:
            logger.print('printing one interacting pair i-j per line\n')
            md.print_verlet_list(self.verlet_list.data)
    
    def zero_accelerations(self,verbose=False):
        """
        Set accelerations to zero
        """
        if not hasattr(self,'ax'):
            self.__init__(self.size(),arrays='a')
        self.ax.fill(0.0)
        self.ay.fill(0.0)
        self.az.fill(0.0)
        if verbose:
            print("ParticleContainer.zero_accelerations - accelerations set to zero.")

    def compute_interactions_verlet(self,zero_accelerations=True,verbose=False,n_iterations=1, verlet_list=None):
        """
        Compute the interactions as specified by self.verlet_list.
        returns the cputime taken in the fortran subprogram
        """
        if zero_accelerations:
            self.zero_accelerations()
        
        if verlet_list:
            vll = verlet_list
        else:
            vll = self.verlet_list
            
        if verbose:
            Log.Print("ParticleContainer.compute_interactions_verlet - computing interactions ...")
        cputime = md.compute_interactions_verlet( self.rx, self.ry, self.rz
                                                , self.ax, self.ay, self.az, vll.data
                                                , n_iterations
                                                )
        if verbose:
            Log.Print("{} x ParticleContainer.compute_interactions_verlet done in {}s".format(n_iterations,cputime))
        return cputime

#     def compute_interactions_hilbert(self,zero_accelerations=True,verbose=False,n_iterations=1):
#         """
#         Compute the interactions as specified by self.verlet_list.
#         returns the cputime taken in the fortran subprogram
#         """
#         if zero_accelerations:
#             self.zero_accelerations()
#             
#         if verbose:
#             print("ParticleContainer.compute_interactions_verlet - computing interactions ...")
#         cputime = md.compute_interactions_hilbert( self.rx, self.ry, self.rz
#                                                  , self.ax, self.ay, self.az, self.hilbert_list
#                                                  , n_iterations
#                                                  )
#         if verbose:
#             print("{} x ParticleContainer.compute_interactions_verlet done in {}s".format(n_iterations,cputime))
#         return cputime
    
#     def r2(self,i,j):
#         ijk_i = np.empty((3,),dtype=np.int32)
#         h=int(self.h[i])
#         print(type(h))
#         print(type(ijk_i))
#         hilbert.h2ijk_1(h,ijk_i)
#         print(i,self.h[i],ijk_i,)
#         print(self.rx[i],self.ry[i],self.rz[i])
# 
#         ijk_j = np.zeros((3,),dtype=np.int32)
#         h=int(self.h[j])
#         hilbert.h2ijk_1(h,ijk_j)
#         print(j,self.h[j],ijk_j)
#         print(self.rx[j],self.ry[j],self.rz[j])
#         r2 = (self.rx[i]-self.rx[j])**2 + (self.ry[i]-self.ry[j])**2 + (self.rz[i]-self.rz[j])**2
#         print(r2)
#         
#         return r2
#--------------------------------------------------------------------------------------------------
# test code below
#--------------------------------------------------------------------------------------------------
    
def testfun_verify_verlet_list(atoms, reserve=50, cell_width=3.0, rcutoff=3.0, verbose=True, spatial_sort=True):
    """
    check that the verlet list built by atoms.build_tables corresponds 
    to the brute force verlet list
    """ 
    assert cell_width>=rcutoff
    if spatial_sort: tf = 'True'
    else:            tf = 'False'
    with Log("testfun_verify_verlet_list(cell_width={},rcutoff={},spatial_sort={})".format(cell_width,rcutoff,tf)) as logger:
        # we must first build the hilbert list because it sorts the atoms
        verlet_list_hilbert = \
        atoms.build_verlet_list('hilbert'    ,cell_width=rcutoff,reserve=reserve,rcutoff=rcutoff,spatial_sort=spatial_sort)
        
        atoms.build_verlet_list('brute_force'                   ,reserve=reserve,rcutoff=rcutoff)
        
#         verlet_list_hilbert.print2(atoms.verlet_list)
        
        missing = verlet_list_hilbert.compare(atoms.verlet_list)
        if missing:
            rcutoff2 = rcutoff*rcutoff
            for m in range(len(missing)):
                if m==0:
                    s = atoms.verlet_list.algo
                else:
                    s = verlet_list_hilbert.algo
                    
                if missing[m]:
                    ok = False
                    for ij in missing[m]:
                        i = ij[0]
                        j = ij[1]
                        r2 = (atoms.rx[i]-atoms.rx[j])**2 + (atoms.ry[i]-atoms.ry[j])**2 + (atoms.rz[i]-atoms.rz[j])**2
                        miss = '!!! missing interaction {}-{} in {}: r2={} rc2={}'.format(i,j,s,r2,rcutoff2)
                        logger.print(miss)
            assert False, 'There are missing interaction pairs in the verlet lists.'
        
def testfun_test_interactions_verlet(atoms):
    """
    Test the loop over interacting pairs using the verlet list
    """
    with Log("testfun_test_interactions_verlet",progressBar=ProgressBar()) as logger:
        atoms.__init__(atoms.size(),arrays='a')
        md.test_interactions_verlet \
            ( atoms.rx, atoms.ry, atoms.rz
            , atoms.ax, atoms.ay, atoms.az, atoms.verlet_list.data
            )
        n_atoms = atoms.size()
    #         print(atoms.ax)
    #         print(atoms.ay)
    #         print(atoms.az)
        for i in range(n_atoms):
            n_pairs_i = atoms.verlet_list.data[0,i]
            n_i_in_other_lists = 0
            for j in range(n_atoms):
                if j!=i:
                    if atoms.verlet_list.find1(j,i):
                        n_i_in_other_lists += 1
    #             print(i,n_pairs_i,n_i_in_other_lists)
            assert atoms.ax[i]==n_pairs_i
            assert atoms.ay[i]==n_pairs_i+n_i_in_other_lists
            assert atoms.az[i]==n_pairs_i-n_i_in_other_lists
    #             print(i,atoms.ay[i],atoms.az[i])

if __name__=="__main__":
    from fcc import FCC
    import math

    ##############################
    # 8 unit cells
    ##############################
    with Log("TEST 8 unit cells"):
        atoms = ParticleContainer()
        assert atoms.size()==0
        
        n_atoms = 4*2**3
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,offset=[0.25,0.25,0.25],verbose=True)
        assert atoms.size()==n_atoms
        assert nb==2
        
        atoms.compute_ijk_h()
        
        i = atoms.i
        j = atoms.j
        k = atoms.k
        del atoms.i
        del atoms.j
        del atoms.k
    
        atoms.compute_ijk()
        assert np.all(i==atoms.i)
        assert np.all(j==atoms.j)
        assert np.all(k==atoms.k)
        
        testfun_verify_verlet_list(atoms,spatial_sort=False)
        testfun_verify_verlet_list(atoms,spatial_sort=False,cell_width=1,rcutoff=1)
        testfun_verify_verlet_list(atoms)
        atoms.verlet_list.print()
        assert atoms.verlet_list.find2(0,1)==1
        assert atoms.verlet_list.find2(0,2)==2
        assert atoms.verlet_list.find2(1,2)==2
        assert atoms.verlet_list.find2(0,3)==3
        assert atoms.verlet_list.find2(1,3)==3
        assert atoms.verlet_list.find2(2,3)==3
        for i in range(n_atoms):
            for j in range(i):
                assert not atoms.verlet_list.find2(i,j) is None, 'Oops: {}-{} not found'.format(i,j)
        assert atoms.verlet_list.n_pairs()==n_atoms*(n_atoms-1)//2
        
        testfun_verify_verlet_list(atoms,reserve=12,rcutoff=math.sqrt(2)/2)
    
        atoms.__init__(n_atoms,arrays='a')
        atoms.print_verlet_list_in_fortran_()
        testfun_test_interactions_verlet(atoms)

    
    ##############################
    # 27 unit cells
    ##############################
    with Log("TEST 27 unit cells"):
        atoms = ParticleContainer()
        assert atoms.size()==0
         
        n = 4*3**3
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n,offset=[0.25,0.25,0.25],verbose=True)
        assert atoms.size()==n
        assert nb==3
         
#         atoms.compute_ijk_h()
        rcutoff = math.sqrt(2)/2
        cell_width =rcutoff
        testfun_verify_verlet_list(atoms,spatial_sort=False,reserve=108)
        testfun_verify_verlet_list(atoms,spatial_sort=False,reserve=108,cell_width=1,rcutoff=1)
        testfun_verify_verlet_list(atoms,reserve=108,cell_width=1,rcutoff=1)
#         testfun_verify_verlet_list(atoms,reserve=12,cell_width=cell_width,rcutoff=rcutoff)
#         testfun_test_interactions_verlet(atoms)
        
    ##############################
    # 1 unit cell
    ##############################
    with Log("TEST 1 unit cell"):
        atoms = ParticleContainer()
        assert atoms.size()==0
        n_atoms = 4
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,offset=[0.25,0.25,0.25],verbose=True)
        assert atoms.size()==4
        assert nb==1
         
        atoms.compute_ijk_h()
        assert np.all(atoms.i==0)
        assert np.all(atoms.j==0)
        assert np.all(atoms.k==0)
        assert np.all(atoms.h==0)
         
        testfun_verify_verlet_list(atoms,spatial_sort=False)
        testfun_verify_verlet_list(atoms)
        atoms.verlet_list.print()
        assert atoms.verlet_list.find2(0,1)==1
        assert atoms.verlet_list.find2(0,2)==2
        assert atoms.verlet_list.find2(1,2)==2
        assert atoms.verlet_list.find2(0,3)==3
        assert atoms.verlet_list.find2(1,3)==3
        assert atoms.verlet_list.find2(2,3)==3
        assert atoms.verlet_list.n_pairs()==6
        atoms.__init__(n_atoms,arrays='a')
        atoms.print_verlet_list_in_fortran_()
        testfun_test_interactions_verlet(atoms)
    
    ##############################
    # 1 unit cell
    ##############################
    with Log("TEST2 1 unit cell"):
        atoms = ParticleContainer()
        assert atoms.size()==0
        n_atoms = 4
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,offset=[0.25,0.25,0.25],verbose=True)
        assert atoms.size()==4
        assert nb==1
        
        atoms.build_verlet_list(cell_width=1,rcutoff=1,verbose=True)
        atoms.verlet_list.print()
        assert np.all(atoms.h==0)

    ##############################
    # 8 unit cells
    ##############################
    with Log("TEST2 8 unit cells"):
        atoms = ParticleContainer()
        assert atoms.size()==0
        
        n_atoms = 4*2**3
        atoms.rx,atoms.ry,atoms.rz,nb = FCC(n_atoms,offset=[0.25,0.25,0.25],verbose=True)
        assert atoms.size()==n_atoms
        assert nb==2
    
        atoms.build_verlet_list(cell_width=1,rcutoff=1,verbose=True)
        for i in range(1,n_atoms):
            assert atoms.h[i]>=atoms.h[i-1]
    #     print(atoms.h)
        h = atoms.h.copy()
        atoms.compute_h()
        assert np.all(h==atoms.h)
    #     print(atoms.hilbert_list)
        for h in range(n_atoms//4):
            assert atoms.hilbert_list[h,0] == 4*h
            assert atoms.hilbert_list[h,1] == 4

    print('done')
    