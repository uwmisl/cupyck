import pandas as pd
import nupyck
from multiprocessing import Pool, cpu_count
from session import Session


def _nupyck_pfunc(kwargs):
    return nupyck.pfunc(**kwargs)


class MulticoreSession(Session):

    def __init__(self, **kwargs):
        ncores = cpu_count()
        options = {
            "ncores": ncores
        }
        if options['ncores'] <= 0 or options['ncores'] > ncores:
            raise ValueError("ncores must be > 0 and <= %d" % ncores)

        super(MulticoreSession, self).__init__(**options)

        self.pool = Pool(self.options['ncores'])

    def shutdown(self):
        self.pool.close()

    def pfunc(self, jobs):

        options = nupyck.Options(
            material = self.options['material'],
            na = self.options['na'],
            mg = self.options['mg'],
            dangles = self.options['dangle_type']
        )

        pfunc_arg_sets = [
            { "sequences": job.sequences,
              "permutation": job.permutation,
              "temp": job.temperature,
              "options": options
            }
            for job in jobs.itertuples()
        ]

        results = pd.DataFrame(
            self.pool.map(_nupyck_pfunc, pfunc_arg_sets)
        )
        results.rename(columns={"pfunc":"partition_function"}, inplace=True)

        return jobs.join(results)
