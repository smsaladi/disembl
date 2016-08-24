"""
This file represents the executeable version of disembl.c wrapped in pure
 python. Use of this as a module requires that both disembl.h and disembl.o
 be in the same directory as this wrapper

Nadine Bradbury - bradbury@caltech.edu

"""

import numpy as np
import ctypes as c
import numpy.ctypeslib as npct

# find and load the disembl.so executeable file

cfile = c.cdll.LibraryLoader('libdisembl.so')
cfile.argstype = ['c_char_p', 'array_1d_float', 'array_1d_float',
                'array_1d_float']

def predict_seq(seq, sm, sb, sr):
    # make a ctypes string of seq
    cseq = c.c_char_p(seq)

    # create the array of pointers to feed into the c function
    sm_arr = npct.ndpointer(dtype=np.float, ndim=1, flags='CONTIGUOUS')
    sb_arr = npct.ndpointer(dtype=np.floar, ndim=1, flags='CONTIGUOUS')
    sr_arr = npct.ndpointer(dtype=np.float, ndim=1, flags='CONTIGUOUS')

    cfile(seq, sm_arr, sb_arr, ar_arr)

    #  pointer arrays back into python lists for calc_disembl
    sm_float = ctypes.cast(sm_arr, ctypes.POINTER(ctypes.c_float))
    sm_list = [sm_float[i] for i in range(arrayLength)]

    sb_float = ctypes.cast(sb_arr, ctypes.POINTER(ctypes.c_float))
    sb_list = [sb_float[i] for i in range(arrayLength)]

    sr_float = ctypes.cast(sr_arr, ctypes.POINTER(ctypes.c_float))
    sr_list = [sr_float[i] for i in range(arrayLength)]

    return (sm_list, sb_list, sr_list)
