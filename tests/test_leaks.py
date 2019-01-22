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

sample_results = []
for i in range(5):
    print "iteration %d" % i
    print "initializing..."

    with cupyck.GPUSession(max_seqlen = 100) as sess:
        print "done"

        print "running pfunc test..."
        results = sess.pfunc(jobs)
        sample_results.append(results.loc[100])
        results = sess.pfunc(jobs)
        sample_results.append(results.loc[100])
        print "done"


print pd.concat(sample_results)['energy']
