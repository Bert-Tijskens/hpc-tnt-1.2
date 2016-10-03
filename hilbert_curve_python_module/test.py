import numpy as np
import imp
import pyHilbert as hilbert
imp.reload(hilbert)
# print(hilbert.__dict__)

if __name__=="__main__":
    print(hilbert.info())
    maxijk= hilbert.cell_index_limit()
    maxh  = hilbert.hilbert_index_limit() 
    print(maxijk)
    print(maxh)
    
    try:
        hilbert.validate_cell_index(maxijk)
    except RuntimeError as e:
        print(e)
    try:
        hilbert.validate_hilbert_index(maxh)
    except RuntimeError as e:
        print(e)
    
    
    for BITS in range(10):
        n = 2**BITS 
        print("BITS =",BITS)
        print("i/j/k <",n)
        print("    h <",n**3)
    
    n=hilbert.cell_index_limit()
    n3=n**3
    print('n =',n)
    print('n3=',n3)
    hh=-np.ones((n3,),dtype=np.int32)
    ijk = np.empty((3,),dtype=np.int32)
    
    for i in range(n):
        for j in range(n):
            for k in range(n):
                #compute h
                h = hilbert.ijk2h_1(i,j,k)
                if not hilbert.is_validating():
                    assert h<n3
                assert hh[h]==-1
                hh[h]=h
                #and from h back to ijk
                ijk.fill(-1) 
                hilbert.h2ijk_1(h,ijk)
                if not hilbert.is_validating():
                    assert i==ijk[0]
                    assert j==ijk[1]
                    assert k==ijk[2]
    print('n =',n)
    print('n3=',n3)
    print("finished.")