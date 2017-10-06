import ctypes
import os
import math
import numpy as np

DNA = 0
RNA = 1

class Session(object):

    lib = None

    def __init__(
            self,
            nblocks,
            nthreads,
            max_seqlen,
            t_lo, t_hi, t_step,
            na=1.0, mg=0.0,
            long_helix=0,
            dangle_type=1,
            material=DNA
            ):

        if t_hi < t_lo:
            raise ValueError("hi temp cannot be less than lo temp")

        if t_step < 0:
            raise ValueError("temp step must be positive")

        if na < 0 or mg < 0:
            raise ValueError("salt concentrations must be positive")

        package_dir = os.path.dirname(__file__)
        lib_dir = os.path.join(package_dir, "../lib")
        os.environ["NUPACKHOME"] = lib_dir

        self.lib = ctypes.cdll.LoadLibrary(
            os.path.join(lib_dir, "pfunc-cuda.so")
        )

        self.lib.pfuncInitialize(
                ctypes.c_int(nblocks),
                ctypes.c_int(nthreads),
                ctypes.c_int(max_seqlen),
                ctypes.c_double(t_lo),
                ctypes.c_double(t_hi),
                ctypes.c_double(t_step),
                ctypes.c_double(na),
                ctypes.c_double(mg),
                ctypes.c_int(long_helix),
                ctypes.c_int(dangle_type),
                ctypes.c_int(material)
                )


    def pfunc(self, seqs, temps, symmetries):
        kB = 0.0019872041
        nseqs = len(seqs)

        seqs = (ctypes.c_char_p * nseqs)(
            *[ ctypes.c_char_p(seq) for seq in seqs ]
        )

        temps = (ctypes.c_double * nseqs)(*temps)

        symmetries = (ctypes.c_int * nseqs)(*symmetries)

        pfs = (ctypes.c_double * nseqs)()

        self.lib.pfuncMulti(
            seqs,
            ctypes.c_int(nseqs),
            symmetries,
            temps,
            ctypes.byref(pfs)
        )

        energies = [
            -kB * (273.15 + temp) * math.log(max(pf, 1))
            for pf, temp in zip(pfs, temps)
        ]

        pfs = np.array(pfs)
        energies = np.array(energies)

        return { "pfs" : pfs, "energies" : energies }
