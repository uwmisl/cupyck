import abc
import itertools
import nupyck.apps.concentrations as concs
import pandas as pd
import numpy as np


DNA = 0
RNA = 1


class Session(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, **kwargs):

        options = {
            "na"         : 1.0,
            "mg"         : 0.0,
            "long_helix" : 1,
            "dangle_type": 1,
            "material"   : DNA
        }
        options.update(kwargs)

        if options['na'] < 0 or options['mg'] < 0:
            raise ValueError("salt concentrations must be non-negative")

        if options['long_helix'] not in [0,1]:
            raise ValueError("invalid long helix specification")

        if options['dangle_type'] not in [0,1,2]:
            raise ValueError("invalid dangle type specification")

        if options['material'] not in [RNA, DNA]:
            raise ValueError("invalid material specification")

        self.options = options

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.shutdown()

    @abc.abstractmethod
    def shutdown(self):
        pass

    @abc.abstractmethod
    def pfunc(self, jobs):
        pass

    def concentrations(self, jobs):

        jobs_with_ids = jobs.join(
            pd.Series(range(len(jobs)), name="job_id", index=jobs.index)
        )

        perms = []
        for job in jobs_with_ids.itertuples():
            for n in range(1, job.max_complex_size + 1):
                job_perms = concs.makePermutations(n, len(job.sequences))

                for perm in job_perms:
                    perms.append((job.job_id, perm))

        perms = pd.DataFrame(perms, columns=["job_id", "permutation"])

        pf_jobs = pd.merge(jobs_with_ids, perms, on="job_id")

        pf_results = self.pfunc(pf_jobs)

        conc_results = []
        for job_id in range(len(jobs)):

            job = jobs.iloc[job_id]
            pf_result = pf_results[pf_results.job_id == job_id]
            perms = pf_result.permutation

            x0   = np.array(job.x0)
            G    = pf_result.energy
            A    = concs._convert_perms_to_A(perms)
            temp = job.temperature

            x = concs.calc_conc(x0, G, A, temp)

            conc_results.append({
                "concentrations": dict(zip(perms, x)),
                "energies": dict(zip(perms, G))
            })

        results = jobs.join(pd.DataFrame(conc_results, index=jobs.index))

        return results

