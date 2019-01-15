import context
import cupyck

import pandas as pd
import numpy as np

randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

jobs = pd.DataFrame([
    { "sequences": [ randseq(20), randseq(20) ],
      "permutation": [1,2],
      "temperature": np.random.choice(range(20,50))
    }
    for _ in range(1000)
])

print "initializing..."
sess = cupyck.GPUSession(max_seqlen = 50)
print "done"

print "running pfunc test..."
results = sess.pfunc(jobs)
print "done"

print "example output:"
print results.loc[0]
