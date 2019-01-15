import context
import cupyck
import multiprocessing

import numpy as np
import pandas as pd

randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

class PFServer(cupyck.Server):

    def worker(self, jobs):
        return self.session.pfunc(jobs)

def server():
    session = cupyck.GPUSession(max_seqlen = 40, nblocks = 100)
    server = PFServer(port = 2046, session = session)
    server.listen(verbose = True)

p = multiprocessing.Process(target = server)
p.start()

jobs = pd.DataFrame(
    { "sequences": [randseq(10), randseq(10)],
      "temperature": np.random.choice(range(20,40)),
      "permutation": [1,2]
    }
    for _ in range(100)
)

calc_pf = cupyck.Client([("localhost", 2046)])

results = calc_pf(jobs)

print "example output:"
print results.loc[0]

p.terminate()
