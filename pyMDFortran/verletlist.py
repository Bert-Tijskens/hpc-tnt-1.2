import numpy as np
from logger import Log,ProgressBar

import sys
sys.path.append('.')

from pyMDFortran import md as md

class VerletList:
    def __init__(self,n_atoms, algo='',reserve=50,rcutoff=None):
        self.data = np.zeros( (reserve+1,n_atoms), dtype=np.int32, order='F')
        self.linear = None
        # mind the fortran storage order! So, we have one column per atom.
        self.algo = algo
        self.rcutoff = rcutoff
    
    def add(self,i_atom,j_atom,test_duplicates=True):
        if test_duplicates:
            if not self.find2(i_atom, j_atom) is None:
                print('Not added (duplicate) {}-{}'.format(i_atom, j_atom))
                return
        if j_atom<i_atom:
            self.data[0                  ,i_atom] += 1
            self.data[self.data[0,i_atom],i_atom] = j_atom
#             print(j_atom,i_atom)
        else:
            self.data[0                  ,j_atom] += 1
            self.data[self.data[0,j_atom],j_atom] = i_atom
#             print(i_atom,j_atom)

    def find2(self,i_atom,j_atom):
        """
        find if pair i_atom,j_atom in somewher in the verlet_lists
        returns
            i_atom if j_atom is in i_atom's verlet list
            j_atom if i_atom is in j_atom's verlet list
            None otherwise
        """
        if self.find1(i_atom,j_atom):
            return i_atom
        elif self.find1(j_atom,i_atom):
            return j_atom
        else:
            return None
        
    def find1(self,i_atom,j_atom):
        """
        return True/False if j_atom is/is not in i_atom's Verlet list
        performs linear search
        """
        n = self.data[0,i_atom]
        for i in range(1,n+1):
            if self.data[i,i_atom]==j_atom:
                return True
        else:
            return False
        
    def n_pairs(self):
        """
        return the total number of entries in the verlet list
        """
        return np.sum(self.data[0,:])
    
    def shape(self):
        return self.data.shape

    def n_atoms(self):
        return self.data.shape[1]

    def output(self):
        """
        Print the verlet list 
        """
        Log.Output("\nVerletList object: algo={}, rcutoff={}".format(self.algo,self.rcutoff))
        for a in range(self.n_atoms()):
            na = self.data[0,a]
            Log.Output('{} {}'.format(a,self.data[0:na+1,a]))
    
    def compare(self,other,reverse=True):
        """
        compare this verlet list with another, 
        return None if both verlet lists are equivalent, and lists with missing interactions otherwise.
        This method may call if reverse is True.
        """
        # self-other
        if reverse:
            message = "VerletList.compare(self,other)"
        else:
            message = "VerletList.compare(other,self) {reversed order}"
        progress = 2+self.n_atoms()//100
        missing = []
        with Log(message,progressBar=ProgressBar()) as logger:
            logger.output("self .VerletList.compare: {} has {} pairs.".format( self.algo, self.n_pairs()))
            logger.output("other.VerletList.compare: {} has {} pairs.".format(other.algo,other.n_pairs()))
            n_atoms = self.n_atoms()
            for i in range(n_atoms):
                logger.show_progress(i/n_atoms)
                n_pairs = self.data[0,i]
                for ji in range(1,n_pairs+1):
                    j = self.data[ji,i]
                    found = other.find2(i,j) 
                    if found is None:
                        missing.append([i,j])
                        s = '!!! missing interaction {}-{}'.format(i,j)
                        logger.output(s)
                        logger.output_final(s)
#                     logger.output('interaction {}-{}'.format(i,j))
                                                
            logger.output_final('n_errors = '+str(len(missing)))
        
        if reverse:
            n_pairs_different = (self.n_pairs()!=other.n_pairs())
            if missing or n_pairs_different:
                # other-self
                return [missing,other.compare(self,reverse=False)]
            else:
                if missing:
                    return [missing,[]]
                else:
                    return None
        else:
            return missing;
        
    def linearize(self):
        """
        transform into linear data structure, removing the unused entries
        [n_pairs_i1 j1_i1 j2_i1 .. n_pairs_i2 j1_i2 j2_i2 .. i3 ..]   
        """
        n_atoms = self.n_atoms()
        n = n_atoms + self.n_pairs()
        self.linear = np.empty((n,),dtype=np.int32)
        j = 0
        for ia in range(n_atoms):
            ia_pairs       = self.data[0,ia]
            self.linear[j] = ia_pairs 
            if ia_pairs:
                self.linear[j+1:j+1+ia_pairs] = self.data[1:1+ia_pairs,ia]
            j += 1+ia_pairs

    def print2(self,other):
        """
        Print 2 verlet lists side by side
        """ 
        print("\nVerletList object: algo=", self.algo)
        print("\nVerletList object: algo=",other.algo)
        for a in range(self.n_atoms()):
            na = self.data[0,a]
            print(a,self.data[0:na+1,a])
            nb = other.data[0,a]
            print(a,other.data[0:nb+1,a])
        
#--------------------------------------------------------------------------------------------------
# test code below
#--------------------------------------------------------------------------------------------------
def testfun_test_interactions_verlet(vl):
    """
    Test the loop over interacting pairs using the verlet list
    """
    import math
    with Log("verletlist.testfun_test_interactions_verlet",progressBar=ProgressBar()) as logger:
        n = vl.n_atoms()
        
        rx = np.zeros((n,),dtype=np.float64)
        ry = np.zeros((n,),dtype=np.float64)
        rz = np.zeros((n,),dtype=np.float64)
        ax = np.zeros((n,),dtype=np.float64)
        ay = np.zeros((n,),dtype=np.float64)
        az = np.zeros((n,),dtype=np.float64)

        md.test_interactions_verlet( rx, ry, rz, ax, ay, az, vl.data )
        for ia in range(n):
            ia_pairs = vl.data[0,ia]
            ia_pairs_in_other_lists = 0
            for j in range(n):
                if j!=ia:
                    if vl.find1(j,ia):
                        ia_pairs_in_other_lists += 1
    #             print(ia,ia_pairs,ia_pairs_in_other_lists)
            assert ax[ia]==ia_pairs
            assert ay[ia]==ia_pairs+ia_pairs_in_other_lists
            assert az[ia]==ia_pairs-ia_pairs_in_other_lists
    #             print(ia,ay[ia],az[ia])

        ax.fill(0.)
        ay.fill(0.)
        az.fill(0.)
        
        if vl.linear is None:
            vl.linearize()
        md.test_interactions_verlet_linear( rx, ry, rz, ax, ay, az, vl.linear )
        j = 0
        for ia in range(n):
            ia_pairs = vl.linear[j]
            ia_pairs_in_other_lists = 0
            for ja in range(n):
                if ja!=ia:
                    if vl.find1(ja,ia):
                        ia_pairs_in_other_lists += 1
    #             print(ia,ia_pairs,ia_pairs_in_other_lists)
            assert ax[ia]==ia_pairs
            assert ay[ia]==ia_pairs+ia_pairs_in_other_lists
            assert az[ia]==ia_pairs-ia_pairs_in_other_lists
    #             print(i,ay[i],az[i])
            j += 1+ia_pairs

        rx.fill(0.)
        ry.fill(0.)
        rz.fill(0.)
        ax.fill(0.)
        ay.fill(0.)
        az.fill(0.)
        md.compute_interactions_verlet_linear( rx, ry, rz, ax, ay, az, vl.linear, 1 )
        print(ax)
        print(ay)
        print(az)
        #since all atoms are at the origin, interatomic distances are zero and forces are infinite (nan)
        j=0
        for ia in range(n):
            ia_pairs = vl.linear[j]
            ia_pairs_in_other_lists = 0
            if ia_pairs==0:
                for ja in range(n):
                    if ja!=ia:
                        if vl.find1(ja,ia):
                            ia_pairs_in_other_lists += 1
            if ia_pairs or ia_pairs_in_other_lists:
                # atom ia participates in at least one interaction henc axyz[ia] is nan 
                assert math.isnan(ax[ia]) 
                assert math.isnan(ay[ia]) 
                assert math.isnan(az[ia])
            else:
                assert ax[ia]==0 
                assert ay[ia]==0 
                assert az[ia]==0 
            j += 1+ia_pairs

if __name__=='__main__':
    # perform tests
    vl = VerletList(8,algo='self')
    assert vl.n_pairs()==0
    vl.output()
    vl.add(0,1)
    vl.add(0,2)
    vl.add(2,1)
    vl.output()
    assert vl.n_pairs()==3
    assert vl.find2(1,2)==2
    assert vl.find2(2,1)==2
    
    vl_other = VerletList(8,algo='other')
    vl_other.add(0,1)
    vl_other.add(0,2)
    vl_other.add(2,1)
    assert vl.compare(vl_other) is None

    vl_other.add(3,1)
    vl_other.output()
    a = vl.compare(vl_other)
    assert vl.compare(vl_other) == [[],[[3,1]]] 
    
    vl_other.linearize()
    assert np.all(vl_other.linear==[0
                                   ,1,0
                                   ,2,0,1 
                                   ,1,1
                                   ,0
                                   ,0
                                   ,0
                                   ,0
                                   ]);
                                   
    md.print_verlet_list  (vl_other.data)
    md.print_verlet_linear(vl_other.n_atoms(), vl_other.linear)
    
    testfun_test_interactions_verlet(vl_other)
    
    testfun_test_interactions_verlet(vl)

    print('done')
     