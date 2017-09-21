import nupyck.core
import nupyck.apps.pfunc
from nupyck.core import DNA
from multiprocessing import Pool
import math

from ctypes import *
lib = cdll.LoadLibrary("./pfunc-cuda.so")

fp_type = c_float

def init():
    lib.pfuncInitialize(
                c_int(16384),
                fp_type(0), fp_type(100), fp_type(0.05),
                fp_type(1.0),
                fp_type(0.0),
                c_int(0),
                c_int(1),
                c_int(0)
            )

def pfmulti(seqs, temps):
    nseqs = len(seqs)

    inputSeqs = (c_char_p * nseqs)(
                *[ c_char_p(seq) for seq in seqs ]
            )

    permSyms = (c_int * nseqs)(*[1 for _ in seqs])

    c_temps = (fp_type * nseqs)(*temps)

    result = (fp_type * nseqs)()
    lib.pfuncMulti(inputSeqs, c_int(nseqs), permSyms, c_temps, byref(result))
    return [-nupyck.core.kB * (273.15 + temp) * math.log(max(pf,1))
            for pf, temp in zip(result, temps)]

def check((seq, temp)):
    return nupyck.apps.pfunc.single(seq, temp, material = DNA)['energy']

def test_seqs():
    import numpy as np
    randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

    seqs = [ randseq(40) + "+" + randseq(40) for _ in range(4096) ]
    temps = np.random.uniform(0,100,4096)

    pool = Pool()
    checks = pool.map(check, zip(seqs, temps))
    pool.close()

    result = pfmulti(seqs, temps)
    print np.abs((np.array(checks) - np.array(result))).max()

init()
test_seqs()
