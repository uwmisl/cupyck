import ctypes
import os
import math
import subprocess
import xml.etree.ElementTree as et

DNA = 0
RNA = 1

class Options(object):
    def __init__(self, max_seqlen, **kwargs):

        options = {
            "nblocks"    : None,
            "nthreads"   : 64,
            "t_lo"       : 0,
            "t_hi"       : 100,
            "t_step"     : 0.05,
            "na"         : 1.0,
            "mg"         : 0.0,
            "long_helix" : 1,
            "dangle_type": 1,
            "material"   : DNA
        }
        options.update(kwargs)

        self.max_seqlen = max_seqlen
        for opt, val in options.items():
            self.__setattr__(opt, val)

        if self.max_seqlen <= 0:
            raise ValueError("max. sequence length must be positive")

        if self.t_hi < self.t_lo:
            raise ValueError("hi temp cannot be less than lo temp")

        if self.t_step <= 0:
            raise ValueError("temp step must be positive")

        if self.na < 0 or self.mg < 0:
            raise ValueError("salt concentrations must be non-negative")

        if self.long_helix not in [0,1]:
            raise ValueError("invalid long helix specification")

        if self.dangle_type not in [0,1,2]:
            raise ValueError("invalid dangle type specification")

        if self.material not in [RNA, DNA]:
            raise ValueError("invalid material specification")

        if self.nblocks is None:
            gpu_info = et.fromstring(
                subprocess.check_output(
                    ["nvidia-smi", "-i", "0", "-q", "-x"]
                )
            )

            mem_free_text = gpu_info.find(".//fb_memory_usage/free").text
            if "MiB" not in mem_free_text:
                raise ValueError("unrecognized memory utilization type")

            mem_free = int(mem_free_text.split(" ")[0]) * (1024 * 1024)

            mem_needed = (
                (max_seqlen ** 2) / 2 * # entries in each pf array
                8                     * # bytes per entry
                10                      # number of arrays
            )

            self.nblocks = min(mem_free / mem_needed, 65536)


class Session(object):

    lib = None

    def __init__(self, options):

        package_dir = os.path.dirname(__file__)
        lib_dir = os.path.join(package_dir, "../lib")
        os.environ["NUPACKHOME"] = lib_dir

        self.lib = ctypes.cdll.LoadLibrary(
            os.path.join(lib_dir, "pfunc-cuda.so")
        )

        self.lib.pfuncInitialize(
                ctypes.c_int(options.nblocks),
                ctypes.c_int(options.nthreads),
                ctypes.c_int(options.max_seqlen),
                ctypes.c_double(options.t_lo),
                ctypes.c_double(options.t_hi),
                ctypes.c_double(options.t_step),
                ctypes.c_double(options.na),
                ctypes.c_double(options.mg),
                ctypes.c_int(options.long_helix),
                ctypes.c_int(options.dangle_type),
                ctypes.c_int(options.material)
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

        return { "pfs" : pfs, "energies" : energies }
