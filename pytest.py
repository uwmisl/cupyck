import nupyck.core as core
import nupyck.core
import nupyck.apps.pfunc
import nupyck.apps.concentrations as conc
from nupyck.core import DNA
from multiprocessing import Pool
import math
import numpy as np

from ctypes import *
lib = cdll.LoadLibrary("./pfunc-cuda.so")

fp_type = c_float

def init():
    lib.pfuncInitialize(
                c_int(256),
                fp_type(0), fp_type(100), fp_type(0.05),
                fp_type(1.0),
                fp_type(0.0),
                c_int(0),
                c_int(1),
                c_int(0)
            )



def pair_perms((template, primer)):
    return [ template,
             primer,
             "+".join((template, template)),
             "+".join((template, primer)),
             "+".join((primer, primer))
            ]


#def melt_search(templates, primers, tconc, pconc, t_lo, t_hi):
def pfmulti(seqs, temps, syms, normed=False):
    nseqs = len(seqs)

    inputSeqs = (c_char_p * nseqs)(
                *[ c_char_p(seq) for seq in seqs ]
            )

    permSyms = (c_int * nseqs)(*syms)

    c_temps = (fp_type * nseqs)(*temps)

    result = (fp_type * nseqs)()
    lib.pfuncMulti(inputSeqs, c_int(nseqs), permSyms, c_temps, byref(result))
    if normed:
        return np.array([-math.log(max(pf,1))
                for pf, temp in zip(result, temps)])
    else:
        return np.array([-nupyck.core.kB * (273.15 + temp) * math.log(max(pf,1))
                for pf, temp in zip(result, temps)])





def calc_conc(G, t_conc, p_conc, temp, na=1.0, mg=0.0):
  """Return the fraction of template converted to template/primer duplex"""

  kT = (temp + 273.15) * core.kB

  waterDensity = core.nupack.WaterDensity(core.c_double(temp))

  # result array
  x = core.c_double_array([-1,-1,-1,-1,-1])

  A = (core.POINTER(core.c_int) * 2)(
    core.c_array([1, 0, 2, 1, 0]),
    core.c_array([0, 1, 0, 1, 2])
  )
  G            = core.c_double_array(G)
  x0           = core.c_double_array([ t_conc/waterDensity, p_conc/waterDensity ])
  numSS        = core.c_int(2)
  numTotal     = core.c_int(5)
  maxIters     = core.c_int(10000)
  tol          = core.c_double(1e-7)
  deltaBar     = core.c_double(1000)
  eta          = core.c_double(0.125)
  kT           = core.c_double(kT)
  maxNoStep    = core.c_int(50)
  maxTrial     = core.c_int(100000)
  perturbScale = core.c_double(100)
  quiet        = core.c_int(1)
  writeLogFile = core.c_int(0)
  logFile      = core.c_char_p(None)
  h20Density   = core.c_double(waterDensity)
  seed         = core.c_int(0)

  converged = core.nupack.CalcConc(
    x, A, G, x0, numSS, numTotal, maxIters, tol, deltaBar,
    eta, kT, maxNoStep, maxTrial, perturbScale, quiet, writeLogFile,
    logFile, h20Density, seed
  )

  if not converged:
    raise RuntimeError("concentration calculation did not converge")

  # convert back to molarity
  x = np.array(x) * waterDensity

  # final [template-primer duplex] vs initial [template]
  conc_ratio = x[3] / t_conc
  return conc_ratio



def test(n, t_lo, t_hi):
    melts = np.zeros(n)
    rem = np.arange(n)
    converged = np.zeros(n).astype(bool)

    randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

    templates = [ randseq(40) for _ in range(n) ]
    primers = [ randseq(40) for _ in range(n) ]
    print [conc.melting_temp(t, p, 1e-1, 1e-1, 0, 100) for t,p in zip(templates,
            primers)]

    seqs = np.concatenate(map(pair_perms, zip(templates, primers))).reshape(-1,5)
    syms = np.tile([1,1,2,1,2], n).reshape(-1,5)

    # check lo temp
    temps = [t_lo] * (n * 5)
    Gs = pfmulti(seqs.flatten(), temps, syms.flatten(),normed=True).reshape(-1, 5)
    fracs = np.array([ calc_conc(G, 1e-1, 1e-1, t_lo) for G in Gs ])

    mask = (fracs <= 0.5)
    melts[rem[mask]] = t_lo
    converged[rem[mask]] = True
    rem = rem[~mask]
    nrem = len(rem)

    if nrem == 0:
        return melts

    # check hi temp
    temps = [t_hi] * (nrem * 5)
    Gs = pfmulti(seqs[rem].flatten(), temps, syms[rem].flatten(),normed=True).reshape(-1, 5)
    fracs = np.array([ calc_conc(G, 1e-1, 1e-1, t_hi) for G in Gs ])

    mask = (fracs >= 0.5)
    melts[rem[mask]] = t_hi
    converged[rem[mask]] = True
    rem = rem[~mask]
    nrem = len(rem)

    if nrem == 0:
        return melts

    t_hi = np.tile(t_hi, nrem)
    t_lo = np.tile(t_lo, nrem)
    eps = 0.02
    while not np.all(converged):
        temps = t_lo + (t_hi - t_lo) / 2
        Gs = pfmulti(seqs[rem].flatten(), np.repeat(temps,5),syms[rem].flatten(),normed=True).reshape(-1,5)
        fracs = np.array([ calc_conc(G, 1e-1, 1e-1, t) for G,t, in zip(Gs, temps)])

        search_lo = (fracs < 0.5 - eps)
        search_hi = (0.5 + eps < fracs)
        mask = ~(search_lo | search_hi)
        melts[rem[mask]] = temps[mask]

        converged[rem[mask]] = True
        rem = rem[~mask]
        nrem = len(rem)

        t_lo[search_hi] = temps[search_hi]
        t_hi[search_lo] = temps[search_lo]
        t_lo = t_lo[~mask]
        t_hi = t_hi[~mask]


    return melts


init()
melts = test(10, 0.0, 100.0)
print melts


def check((seq, temp)):
    return nupyck.apps.pfunc.single(seq, temp, material = DNA)['energy']

def test_seqs():
    import numpy as np
    randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

    seqs = [ randseq(40) + "+" + randseq(40) for _ in range(4096) ]
    temps = np.random.uniform(0,100,4096)

    pool = Pool()
    checks = pool.map(check, zip(seqs, temps))
    pool.close()

    result = pfmulti(seqs, temps)
    print np.abs((np.array(checks) - np.array(result))).max()

#init()
#test_seqs()
