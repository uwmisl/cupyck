import nupyck.core

seq = "GACAGAGTACGACTAGCATCGATCGACTG"
sint = nupyck.core.c_array(nupyck.core.seqToInts(seq))

from ctypes import *
lib = cdll.LoadLibrary("./pfunc-cuda.so")

init = lib.pfuncInitialize
init(c_double(273.15 + 62), c_double(1.0), c_double(0.0),
        c_int(0), c_int(1), c_int(0))

pfunc = lib.pfuncFullWithSym
pfunc.restype = c_double

pf = pfunc(sint, 1)
print(pf)

import nupyck.apps.pfunc
pf = nupyck.apps.pfunc.single(seq, 62, material=nupyck.apps.pfunc.core.DNA)
print (pf)
