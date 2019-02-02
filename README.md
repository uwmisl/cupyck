cupyck
======
`cupyck` is a python frontend to [NUPACK](nupack.org) that emphasizes data-parallel
execution (running many *independent* NUPACK jobs in parallel). At the core of `cupyck` is
a CUDA-accelerated implementation of NUPACK's partition function, which provides a significant
speedup over the CPU implementation, especially for batch execution.

In addition to frontends for the `pfunc` and `concentrations` applications, `cupyck`'s interface
includes basic support for running jobs on a cluster of remote machines (such as AWS instances).

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

#### Creating a Session

#### `pfunc`

#### `concentrations`

### Advanced Usage

#### Writing Applications

#### Remote Execution
