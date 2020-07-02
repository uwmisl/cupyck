cupyck
======
`cupyck` is a python frontend to [NUPACK](nupack.org) that emphasizes data-parallel
execution (running many *independent* NUPACK jobs in parallel). At the core of `cupyck` is
a CUDA-accelerated implementation of NUPACK's partition function, which provides a significant
speedup over the CPU implementation, especially for batch execution. There is
also basic support for remote execution (such as on a cluster of GPU machines).

`cupyck` has limitations. At the moment, there are only frontends for the
`pfunc` and `concentrations` applications from NUPACK (no `pairs`, `mfe`, etc.),
and computation of RNA pseudoknots is not supported.

Due to issues with floating-point data types, you may find that `cupyck`
produces slightly different numerical results than NUPACK. If this is an issue
for your use case, we recommend that you NUPACK instead.

Table of Contents:
1. [Installation](#installation)
2. [Basic Usage](#basic-usage)
    1. [Creating a Session](#creating-a-session)
    2. [`pfunc`](#pfunc)
    3. [`concentrations`](#concentrations)
3. [Advanced Usage](#advanced-usage)
    1. [Writing Applications](#writing-applications)
    2. [Remote Execution](#remote-execution)

### Installation
Install with `pip install -r requirements.txt .`. This will download the required dependencies, and
compile and install the library. If you wish to use `cupyck` on a GPU-enabled server or workstation,
you must have CUDA and the NVIDIA CUDA compiler (`nvcc`) installed.

### Basic Usage

**Note**: This example is for illustration only. In reality you should be
running `pfunc` with thousands of sequences to amortize the startup costs.
```python
import cupyck
import pandas as pd

jobs = pd.DataFrame(
    [ { "sequences": ["GATCGATCAGAT", "TAGACACATCAG", "TAGAGCAGTA"],
        "permutation": [1,2,2,3],
        "temperature": 23
      },
      { "sequences": ["GGATACAGAT", "TAGAGACCCCG"],
        "permutation": [1,2],
        "temperature": 40
      }
    ]
)

with cupyck.GPUSession(max_seqlen=100) as sess:
    results = sess.pfunc(jobs)

print results.to_html()
```
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>permutation</th>
      <th>sequences</th>
      <th>temperature</th>
      <th>energy</th>
      <th>partition_function</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>[1, 2, 2, 3]</td>
      <td>[GATCGATCAGAT, TAGACACATCAG, TAGAGCAGTA]</td>
      <td>23</td>
      <td>-17.129853</td>
      <td>4.375900e+12</td>
    </tr>
    <tr>
      <th>1</th>
      <td>[1, 2]</td>
      <td>[GGATACAGAT, TAGAGACCCCG]</td>
      <td>40</td>
      <td>-5.220261</td>
      <td>4.397322e+03</td>
    </tr>
  </tbody>
</table>

#### Creating a Session

A "session" holds the context that each call to `pfunc` will execute under (for
instance, GPU memory that only needs to be allocated once). When you create
a session, you specify the hardware parameters (number of GPU or CPU threads,
etc) as well as certain parameters (such as salt concentration) that cannot
change between calls to `pfunc` due to technical limitations.

There are currently two types of sessions, `GPUSession` (for executing parallel
jobs on a CUDA-enabled GPU) and `MulticoreSession` (for executing parallel jobs
on a multicore CPU). Both types of sessions require you to fix the following
parameters at creation time:

| Parameter    | Default Value | Description                       |
| ---------    | ------------- | --------------------------------- |
| `na`         | 1.0           | sodium concentration              |
| `mg`         | 0.0           | magnesium concentration           |
| `long_helix` | 1             | nupack's "long helix" correction" |
| `dangle_type`| 1             | nupack's dangle setting           |
| `material`   | DNA           | energy model (DNA or RNA)         |

##### GPU Sessions
When creating a GPU session, you must specify the following additional
parameters:

| Parameter    | Default Value | Description                       |
| ---------    | ------------- | --------------------------------- |
| `max_seqlen` | must specify  | length (nt) of longest complex    |
| `nblocks`    | inferred      | number of parallel jobs           |
| `nthreads`   | 64            | number of threads per job         |
| `t_lo`       | 0             | lowest temperature to process     |
| `t_hi`       | 100           | highest temperature to process    |
| `t_step`     | 0.05          | precision of temps to process     |

The value of `max_seqlen` must account for the length of
any ordered complex you will be processing. For instance, if you have three
single-stranded species that are 40-nt long, and will be calculating the energy
of a triple-stranded complex, you must specify the maximum sequence length as
120 nucleotides.

If you try to process a temperature between `t_step`s, it will be rounded to the
nearest step.

The `nblocks` parameter, if not specified, is determined by the amount of memory
your GPU hardware has available. `nthreads` is set to a reasonable default, but
should be tuned for your hardware and use case.

##### Multicore Sessions
Multicore sessions only require one additional parameter, `ncores`, which
is the number of CPU cores to use for parallel execution. By default, it is
set to the number of cores the system has.

#### `pfunc`
Once you have created a session, you can call its `pfunc` method on a list of jobs
defined as a [pandas](https://pandas.pydata.org/) data frame. The data frame
must contain the following columns (which mirror the inputs to 
[`nupyck.pfunc`](https://github.com/uwmisl/nupyck#pfunc)):

| Column        | Description                        |
| ------------- | ---------------------------------- |
| `sequences`   | the list of sequences for each job |
| `permutation` | the ordered complex for each job   |
| `temperature` | the temperature for each job       |

The return value is the input data frame augmented with two additional columns:

| Column               | Description                                  |
| -------------------  | -------------------------------------------- |
| `partition_function` | the computed partition function for each job |
| `energy`             | the computed free energy for each job        |

The input data frame may contain additional columns. These will not be used or modified
(unless their name conflicts with one of the output columns), and they will be returned
with the output.

#### `concentrations`

### Advanced Usage

#### Writing Applications

#### Remote Execution
