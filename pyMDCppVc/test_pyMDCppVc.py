from __future__ import division
import pyMDCppVc as md
import numpy as np
import math

def verlet_list_linear_to_nj(vll,n_atoms,N=4):
    """
    convert verlet_list_linear vll to vln and vlj (padded to N)
    N is width of simd registers (better abstraction needed) 
    !!!
        hardcoding N makes that the tests below fail 
        - if the vector register has a different widht,
        - if the pyMDCppVc library was compiled without -mavx (gcc) or -xAVX (icpc) 
    !!!
    """
    vln = np.empty( (n_atoms,), dtype=np.int32 )
    # first compute the size of vlj (which is padded)
    k=0
    n = 0
    for i in range(n_atoms):
        n_pairs_i = vll[k]
        k += n_pairs_i + 1
        while n_pairs_i%N>0: # add an entry for padding
            n_pairs_i += 1
        vln[i] = n_pairs_i
        n+=n_pairs_i
        
#     print(vln)
    vlj = -np.ones( (n,), dtype=np.int32 )
    # next copy the pairs and pad
    k=0
    l=0
    for i in range(n_atoms):
        n_pairs_i        = vll[k]
        n_pairs_i_padded = vln[i]
        k += 1
#         print(n_pairs_i,n_pairs_i_padded)
#         print(l,k)
        if n_pairs_i>0:
            vlj[l:l+n_pairs_i] = vll[k:k+n_pairs_i]
#             print('*',vlj)
            last = vlj[l+n_pairs_i-1]
#             print('**',last)
            vlj[l+n_pairs_i:l+n_pairs_i_padded] = last #copy last entry into padded entries
#             print('***',vlj)
            k += n_pairs_i
            l += n_pairs_i_padded
    return vln,vlj

def compute_interactions_scalar(rx,ry,rz,ax,ay,az,verlet_n,verlet_j,imax=None):
    ax.fill(0)
    ay.fill(0)
    az.fill(0)
    n_atoms = rx.shape[0]
    if imax is None:
        n = n_atoms
    else:
        n = imax+1
    k0 = 0
    for i in range(n):
        n_pairs = verlet_n[i]
        for k in range(n_pairs):
            j = verlet_j[k0+k]
            dx = rx[j]-rx[i]
            dy = ry[j]-ry[i]
            dz = rz[j]-rz[i]
            r2 = dx*dx+dy*dy+dz*dz
            f = ff(r2)
            dx*=f
            dy*=f
            dz*=f
            ax[i]+=dx
            ay[i]+=dy
            az[i]+=dz
            ax[j]-=dx
            ay[j]-=dy
            az[j]-=dz
        m = n_pairs%4
        if m>0:
            m=4-m
        k0 += n_pairs + m

itest = 0
def compute(rx,ry,rz,ax,ay,az,verlet_n,verlet_j,imax=-1):
    global itest
    itest += 1
    print("\ntest {} : input".format(itest))
    print("rx={}".format(rx))
    print("ry={}".format(ry))
    print("rz={}".format(rz))
    ax.fill(0)
    ay.fill(0)
    az.fill(0)
    print("verlet_n={}".format(verlet_n))
    print("verlet_j={}".format(verlet_j))
    
    print("c++: md.compute_interactions_verlet_list(...")
    cputime = md.compute_interactions_verlet_list_test( rx, ry, rz, ax, ay, az, verlet_n, verlet_j, 1, imax ) 
    print("c++: md.compute_interactions_verlet_list(...)")
    print("\ntest {} : output".format(itest))
    print("ax={}".format(ax))
    print("ay={}".format(ay))
    print("az={}".format(az))
    
def alloc(n_atoms,n_pairs):
    rx = np.zeros( (n_atoms,), dtype=np.float64 )
    ry = np.zeros( (n_atoms,), dtype=np.float64 )
    rz = np.zeros( (n_atoms,), dtype=np.float64 )
    ax = np.zeros( (n_atoms,), dtype=np.float64 )
    ay = np.zeros( (n_atoms,), dtype=np.float64 )
    az = np.zeros( (n_atoms,), dtype=np.float64 )
    bx = np.zeros( (n_atoms,), dtype=np.float64 )
    by = np.zeros( (n_atoms,), dtype=np.float64 )
    bz = np.zeros( (n_atoms,), dtype=np.float64 )
    vln= np.zeros( (n_atoms,), dtype=np.int32   )
    if n_pairs>0:
        vlj= np.zeros( (n_pairs,), dtype=np.int32   )
    else:
        vlj = None
    assert(n_pairs%4==0)
    return rx,ry,rz,ax,ay,az,bx,by,bz,vln,vlj

def ff(r2):
        rec_r2 = 1/r2
        rec_r6 = rec_r2*rec_r2*rec_r2
        ff = 48*rec_r2*rec_r6*( rec_r6 - 0.5 )
        return ff;    


if __name__ == "__main__":
    rmin = math.pow( 2, 1/6 )
    print("rmin    ={}".format(rmin))
    print("rmin**2 ={}".format(rmin*rmin))

    # test 1
    rx,ry,rz,ax,ay,az,bx,by,bz,vln,vlj = alloc( 5, 4 ) 
    for i in range(0,4):
        rx[i] = rmin
    vln[4] = 4
    vlj[:] = [0,1,2,3]
    
    compute( rx, ry, rz, ax, ay, az, vln, vlj )
    ffmin = ax[0] 
    assert( np.all( np.abs(ax) < 1e-14 ) )
    assert( ax[0] == ffmin )
    assert( ax[1] == ffmin )
    assert( ax[2] == ffmin )
    assert( ax[3] == ffmin )
    assert( ax[4] == -4*ffmin )
    assert( np.all( np.abs(ay) == 0 ) )
    assert( np.all( np.abs(az) == 0 ) )
    compute_interactions_scalar(rx,ry,rz,bx,by,bz,vln,vlj)
    for i in range(rx.shape[0]):
        assert bx[i]==ax[i]
        assert by[i]==ay[i]
        assert bz[i]==az[i]
    
    # test 2
    vln[4] = 3
    vlj[:] = [0,1,2,2] #last element is dummy 
    compute(rx, ry, rz, ax, ay, az, vln, vlj )
    assert( np.all( np.abs(ax) < 1e-14 ) )
    assert( ax[0] == ffmin )
    assert( ax[1] == ffmin )
    assert( ax[2] == ffmin )
    assert( ax[3] == 0 )
    assert( ax[4] == -3*ffmin )
    assert( np.all( np.abs(ay) == 0 ) )
    assert( np.all( np.abs(az) == 0 ) )
    compute_interactions_scalar(rx,ry,rz,bx,by,bz,vln,vlj)
    for i in range(rx.shape[0]):
        assert bx[i]==ax[i]
        assert by[i]==ay[i]
        assert bz[i]==az[i]

    # test 3
    rx,ry,rz,ax,ay,az,bx,by,bz,vln,vlj = alloc( 6, 8 ) 
    rx[5] = rmin
    vln[5] = 5
    vlj[:] = [0,1,2,3
             ,4,4,4,4] #last 3 element are dummies 
    compute(rx, ry, rz, ax, ay, az, vln, vlj )
    assert( np.all( np.abs(ax) < 1e-13 ) )
    for i in range(5):
        assert( ax[i] ==-ffmin )
    assert ax[5] == 5*ffmin
    assert np.all( np.abs(ay) == 0 ) 
    assert np.all( np.abs(az) == 0 )
    
    compute_interactions_scalar(rx,ry,rz,bx,by,bz,vln,vlj)
    for i in range(rx.shape[0]):
        assert bx[i]==ax[i]
        assert by[i]==ay[i]
        assert bz[i]==az[i]
    
    # test 4
    n_atoms = 10
    rx,ry,rz,ax,ay,az,bx,by,bz,vln,vlj = alloc( n_atoms, 0 )
    for i in range(n_atoms):
        rx[i] = i*rmin
    vlj = []
    for i in range(1,n_atoms):
        vln[i] = i
        for j in range(i):
            vlj.append(j)
        while len(vlj)%4 != 0:
            vlj.append(j)
        print(vlj)
    for imax in range(n_atoms):
        print("imax={}".format(imax))
        compute(rx, ry, rz, ax, ay, az, vln, vlj, imax )
    #     for i in range(n_atoms):
    #         assert(ax[i])
        assert np.all( np.abs(ay) == 0 )
        assert np.all( np.abs(az) == 0 )
    
        compute_interactions_scalar(rx,ry,rz,bx,by,bz,vln,vlj, imax)
        for i in range(imax+1):
            if not bx[i]==ax[i]:
                print("imax={}".format(imax))
                print("i=i{}".format(i))
                print(bx[i])
                print(ax[i])
            assert abs(bx[i]-ax[i])<1e-14
            
            assert by[i]==ay[i]
            assert bz[i]==az[i]
        print("imax={} : ok\n".format(imax))
        
    # test verlet_list_linear_to_nj
    n_atoms = 10
    n_pairs = n_atoms*(n_atoms-1)//2
    vll = -np.ones( (n_atoms + n_pairs,), dtype=np.int32 )
    k=0
    for i in range(n_atoms):
        vll[k] = i 
        k += 1
        for j in range(i):
            vll[k+j] = j 
        k += i
        print(vll)
    k=0
    for i in range(n_atoms):
        n_pairs_i = vll[k]
        print( vll[k:k+1+n_pairs_i])
        k += 1+n_pairs_i
    vln,vlj = verlet_list_linear_to_nj(vll,n_atoms)
    print(vln)
    print(vlj)
    k=0
    for i in range(n_atoms):
        print( i, vln[i], vlj[k:k+vln[i]])
        k += vln[i]
    

    print("test_pyMDCppVc.py done")