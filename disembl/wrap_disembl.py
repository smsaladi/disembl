"""
This file represents the executeable version of disembl.c wrapped in pure
 python. Use of this as a module requires that both disembl.h and disembl.o
 be in the same directory as this wrapper

Nadine Bradbury - bradbury@caltech.edu

"""

import numpy as np
import ctypes

clib = np.ctypeslib.load_library("libdisembl.so", ".")

def pypredict_seq(seq):
    # create arrays for disembl output
    sm_arr = np.zeros(len(seq))
    sb_arr = np.zeros_like(sm_arr)
    sr_arr = np.zeros_like(sm_arr)

    # http://stackoverflow.com/a/37888716/2320823
    clib.predict_seq(ctypes.c_char_p(str.encode(seq)),
                     np.ctypeslib.as_ctypes(sm_arr),
                     np.ctypeslib.as_ctypes(sb_arr),
                     np.ctypeslib.as_ctypes(sr_arr))

    return (sm_arr.tolist(), sb_arr.tolist(), sr_arr.tolist())


if __name__ == "__main__":
    print(pypredict_seq("MAKQPGLDFQSAKGGLGELKRRLLFVIGALIVFRIGSFIPIPGIDAAVLAKLLEQQRGTI"))
