from session import Session
import subprocess
import ctypes
import xml.etree.ElementTree as et
from nupyck.core import calcVPi
import os
import math
import pandas as pd


class GPUSession(Session):

    lib = None

    def __init__(self, max_seqlen, **kwargs):
        try:
            subprocess.check_output("nvidia-smi")
        except:
            raise RuntimeError(
                "No GPU is available."
            )


        options = {
            "max_seqlen" : max_seqlen,
            "nblocks"    : None,
            "nthreads"   : 64,
            "t_lo"       : 0,
            "t_hi"       : 100,
            "t_step"     : 0.05,
        }
        options.update(kwargs)

        if options['nblocks'] is None:
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

            options['nblocks'] = min(mem_free / mem_needed, 65536)

        super(GPUSession, self).__init__(**options)

        backend_dir = os.path.dirname(__file__)
        package_dir = os.path.abspath(os.path.join(backend_dir, '..'))
        os.environ["NUPACKHOME"] = package_dir

        self.lib = ctypes.cdll.LoadLibrary(
            os.path.join(package_dir, "cupyck.so")
        )

        options = self.options
        self.lib.pfuncInitialize(
                ctypes.c_int(options['nblocks']),
                ctypes.c_int(options['nthreads']),
                ctypes.c_int(options['max_seqlen']),
                ctypes.c_double(options['t_lo']),
                ctypes.c_double(options['t_hi']),
                ctypes.c_double(options['t_step']),
                ctypes.c_double(options['na']),
                ctypes.c_double(options['mg']),
                ctypes.c_int(options['long_helix']),
                ctypes.c_int(options['dangle_type']),
                ctypes.c_int(options['material'])
                )

    def shutdown(self):
        pass

    def pfunc(self, jobs):
        kB = 0.0019872041
        njobs = len(jobs)

        seqs = [
            "+".join(
                job.sequences[p-1] for p in job.permutation
            )
            for job in jobs.itertuples()
        ]

        seqs = (ctypes.c_char_p * njobs)(
            *[ ctypes.c_char_p(seq) for seq in seqs ]
        )

        temps = (ctypes.c_double * njobs)(*jobs.temperature)

        symmetries = (ctypes.c_int * njobs)(
            *[ calcVPi(
                   (ctypes.c_int * len(perm))(*perm),
                   ctypes.c_int(len(perm))
               )
               for perm in jobs.permutation
             ]
        )

        pfs = (ctypes.c_double * njobs)()

        self.lib.pfuncMulti(
            seqs,
            ctypes.c_int(njobs),
            symmetries,
            temps,
            ctypes.byref(pfs)
        )

        energies = [
            -kB * (273.15 + temp) * math.log(max(pf, 1))
            for pf, temp in zip(pfs, temps)
        ]

        results = pd.DataFrame(
            zip(energies, pfs),
            columns = ['energy', 'partition_function'],
            index = jobs.index
        )

        return pd.concat([jobs, results], axis=1)

