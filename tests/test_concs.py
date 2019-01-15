import context
import cupyck

import pandas as pd
import numpy as np

randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

jobs = pd.DataFrame([
    { "sequences": [ randseq(20), randseq(20) ],
      "temperature": np.random.choice(range(20,50)),
      "max_complex_size": 2,
      "x0": np.array([1e-6, 1e-6])
    }
    for _ in range(1000)
])

print "initializing..."
sess = cupyck.GPUSession(max_seqlen = 50)
print "done"

print "running concentrations test..."
results = sess.concentrations(jobs)
print "done"

print "example output:"
print results.loc[0]
