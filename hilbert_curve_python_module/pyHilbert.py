# python extensions to the pyHilbertCpp module
from pyHilbertCpp import *
import numpy as np

def reorder(I,told,tnew):
    """
    A python wrapper around the various reorder_<dtype> functions exposed by pyHilbertCpp
    """
    assert isinstance(I   ,np.ndarray),"First argument (I) must be numpy array."
    assert isinstance(told,np.ndarray),"Second argument (t_old) must be numpy array."
    assert isinstance(tnew,np.ndarray),"Third argument (t_new) must be numpy array."
    assert I.shape==told.shape and told.shape==tnew.shape,"Arguments (I, t_old, t_new) must have the same shape."
    assert told.dtype==tnew.dtype,"Arguments t_old and t_new  must have the same dtype."
    
#     print("told: ",told)
    if told.dtype==np.float64:
        print("told.dtype==np.float64")
        reorder_float64(I,told,tnew) 
    elif told.dtype==np.float32:
        print("told.dtype==np.float32")
        reorder_float32(I,told,tnew) 
    elif told.dtype==np.int32:
        print("told.dtype==np.int32")
        reorder_int32(I,told,tnew)
    elif told.dtype==np.uint32:
        print("told.dtype==np.uint32")
        reorder_uint32(I,told,tnew)
    else:
        raise NotImplemented()
#     print("tnew",tnew)