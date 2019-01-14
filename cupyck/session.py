import ctypes
import os
import math
import subprocess
import itertools
import xml.etree.ElementTree as et
import pandas as pd
import nupyck.apps.concentrations as concs

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
        os.environ["NUPACKHOME"] = package_dir

        self.lib = ctypes.cdll.LoadLibrary(
            os.path.join(package_dir, "cupyck.so")
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
            *[ concs.core.calcVPi(
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

    def concentrations(self, jobs):

        jobs_with_ids = jobs.join(pd.Series(range(len(jobs)), name="job_id"))

        perms = []
        for job_id, job in jobs.iterrows():
            for n in range(1, job.max_complex_size + 1):
                job_perms = itertools.combinations_with_replacement(
                    range(1, len(job.sequences) + 1), n
                )
                for perm in job_perms:
                    perms.append((job_id, perm))

        perms = pd.DataFrame(perms, columns=["job_id", "permutation"])

        pf_jobs = pd.merge(jobs_with_ids, perms, on="job_id")

        pf_results = self.pfunc(pf_jobs)

        conc_results = []
        for job_id in range(len(jobs)):

            pf_result = pf_results[pf_results.job_id == job_id]
            perms = pf_result.permutation

            x0   = jobs.loc[job_id].x0
            G    = pf_result.energy
            A    = concs._convert_perms_to_A(perms)
            temp = jobs.loc[job_id].temperature

            x = concs.calc_conc(x0, G, A, temp)

            conc_results.append({
                "concentrations": dict(zip(perms, x)),
                "energies": dict(zip(perms, G))
            })

        results = jobs.join(pd.DataFrame(conc_results))

        return results
