import nupyck.core
import nupyck.apps.pfunc

from ctypes import *
lib = cdll.LoadLibrary("./pfunc-cuda.so")
lib.pfuncInitialize(
            c_int(64),
            c_double(273.15 + 62),
            c_double(1.0),
            c_double(0.0),
            c_int(0),
            c_int(1),
            c_int(0)
        )

def pfmulti(seqs):
    nseqs = len(seqs)

    inputSeqs = (c_char_p * nseqs)(
                *[ c_char_p(seq) for seq in seqs ]
            )

    permSyms = (c_int * nseqs)(*[1 for _ in seqs])

    result = (c_double * nseqs)()
    lib.pfuncMulti(inputSeqs, c_int(nseqs), permSyms, byref(result))
    return result

def test_seqs():
    import numpy as np
    randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

    seqs = [ randseq(40) + "+" + randseq(40) for _ in range(128) ]

    check_pf = lambda seq: nupyck.apps.pfunc.single(seq, 62, material =
            nupyck.apps.pfunc.core.DNA)['pfunc']

    checks = map(check_pf, seqs)

    result = pfmulti(seqs)
    print (np.array(checks) - np.array(result)).max()

test_seqs()
