import numpy as np
from logger import Log,ProgressBar

class VerletList:
    def __init__(self,n_atoms, algo='',reserve=50,rcutoff=None):
        self.data = np.zeros( (reserve+1,n_atoms), dtype=np.int32, order='F')
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
            self.data[0                   ,j_atom] += 1
            self.data[self.data[0,j_atom,],j_atom] = i_atom
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

    def print(self):
        """
        Print the verlet list 
        """
        Log.Print("\nVerletList object: algo={}, rcutoff={}".format(self.algo,self.rcutoff))
        for a in range(self.n_atoms()):
            na = self.data[0,a]
            Log.Print('{} {}'.format(a,self.data[0:na+1,a]))
    
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
            logger.print("self .VerletList.compare: {} has {} pairs.".format( self.algo, self.n_pairs()))
            logger.print("other.VerletList.compare: {} has {} pairs.".format(other.algo,other.n_pairs()))
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
                        logger.print(s)
                        logger.print_final(s)
#                     logger.print('interaction {}-{}'.format(i,j))
                                                
            logger.print_final('n_errors = '+str(len(missing)))
        
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
if __name__=='__main__':
    vl = VerletList(8,algo='self')
    assert vl.n_pairs()==0
    vl.print()
    vl.add(0,1)
    vl.add(0,2)
    vl.add(2,1)
    vl.print()
    assert vl.n_pairs()==3
    assert vl.find2(1,2)==2
    assert vl.find2(2,1)==2
    
    vl_other = VerletList(8,algo='other')
    vl_other.add(0,1)
    vl_other.add(0,2)
    vl_other.add(2,1)
    assert vl.compare(vl_other) is None

    vl_other.add(3,1)
    vl_other.print()
    a = vl.compare(vl_other)
    assert vl.compare(vl_other) == [[],[[3,1]]] 
    
    print('done')
     