import nupyck.core
import nupyck.apps.pfunc
from nupyck.core import DNA
import math

from ctypes import *
lib = cdll.LoadLibrary("./pfunc-cuda.so")
def init():
    lib.pfuncInitialize(
                c_int(64),
                c_double(0), c_double(100), c_double(0.125),
                c_double(1.0),
                c_double(0.0),
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

    c_temps = (c_double * nseqs)(*temps)

    result = (c_double * nseqs)()
    lib.pfuncMulti(inputSeqs, c_int(nseqs), permSyms, c_temps, byref(result))
    return [-nupyck.core.kB * (273.15 + temp) * math.log(max(pf,1))
            for pf, temp in zip(result, temps)]

def test_seqs():
    import numpy as np
    randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

    seqs = [ randseq(80) + "+" + randseq(80) for _ in range(256) ]
    temps = np.random.uniform(0,100,256)

    checks = [ nupyck.apps.pfunc.single(seq, temp, material = DNA)['energy']
            for seq, temp in zip(seqs, temps) ]

    result = pfmulti(seqs, temps)
    print np.abs((np.array(checks) - np.array(result))).max()

init()
test_seqs()
