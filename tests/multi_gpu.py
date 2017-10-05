import zmq
import pickle
from ctypes import *
import numpy as np
import math

from nupyck.apps.concentrations import calc_conc

port = 2046

def pair_perms((template, primer)):
    return [ template,
             primer,
             "+".join((template, template)),
             "+".join((template, primer)),
             "+".join((primer, primer))
            ]

def melting_temps(templates, primers, t_conc, p_conc, t_lo, t_hi):
    t_lo = float(t_lo)
    t_hi = float(t_hi)

    x0 = np.array([t_conc, p_conc])
    A  = np.array([[1, 0, 2, 1, 0], [0, 1, 0, 1, 2]])

    def calc_frac(G, t):
        x = calc_conc(x0, G, A, t)
        frac = x[3] / t_conc
        return frac

    n = len(templates)
    melts = np.zeros(n)
    rem = np.arange(n)
    converged = np.zeros(n).astype(bool)

    seqs = np.concatenate(map(pair_perms, zip(templates, primers))).reshape(-1,5)
    syms = np.tile([1,1,2,1,2], n).reshape(-1,5)

    # check lo temp
    temps = [t_lo] * (n * 5)
    Gs = pfmulti(seqs.flatten(), temps, syms.flatten()).reshape(-1, 5)
    fracs = np.array([ calc_frac(G, t_lo) for G in Gs ])

    mask = (fracs <= 0.5)
    melts[rem[mask]] = t_lo
    converged[rem[mask]] = True
    rem = rem[~mask]
    nrem = len(rem)

    if nrem == 0:
        return melts

    # check hi temp
    temps = [t_hi] * (nrem * 5)
    Gs = pfmulti(seqs[rem].flatten(), temps, syms[rem].flatten()).reshape(-1, 5)
    fracs = np.array([ calc_frac(G, t_hi) for G in Gs ])

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
        Gs = pfmulti(seqs[rem].flatten(), np.repeat(temps,5),syms[rem].flatten()).reshape(-1,5)
        fracs = np.array([ calc_frac(G, t) for G,t, in zip(Gs, temps)])

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


def client(server_list, pairs):

    nservers = len(server_list)
    nseqs = len(pairs)

    splits = range(0, nseqs + 1, nseqs / nservers)
    splits[-1] = nseqs
    ranges = zip(splits, splits[1:])

    chunks = [ pairs[start:end] for start, end in ranges ]

    ctxt = zmq.Context()
    socks = []
    for server, chunk in zip(server_list, chunks):
        sock = ctxt.socket(zmq.REQ)
        sock.connect("tcp://%s:%d" % (server, port))
        sock.send(pickle.dumps(chunk))
        socks.append(sock)

    results = []
    for sock in socks:
        msg = sock.recv()
        result = pickle.loads(msg)
        results.append(result)

    return np.concatenate(results)



lib = None
fp_type = c_double
def init():
    global lib
    lib = cdll.LoadLibrary("../lib/pfunc-cuda.so")

    lib.pfuncInitialize(
                c_int(16384),
                fp_type(0), fp_type(100), fp_type(0.05),
                fp_type(1.0),
                fp_type(0.0),
                c_int(0),
                c_int(1),
                c_int(0)
            )


def pfmulti(seqs, temps, syms):
    kB = 0.0019872041
    nseqs = len(seqs)

    inputSeqs = (c_char_p * nseqs)(
                *[ c_char_p(seq) for seq in seqs ]
            )

    permSyms = (c_int * nseqs)(*syms)

    c_temps = (fp_type * nseqs)(*temps)

    result = (fp_type * nseqs)()
    lib.pfuncMulti(inputSeqs, c_int(nseqs), permSyms, c_temps, byref(result))

    return np.array([-kB * (273.15 + temp) * math.log(max(pf,1))
            for pf, temp in zip(result, temps)])


def server():

    ctxt = zmq.Context()
    sock = ctxt.socket(zmq.REP)
    sock.bind("tcp://*:%d" % port)

    print "initializing..."
    init()

    print "listening."
    while True:
        msg = sock.recv()
        pairs = pickle.loads(msg)
        templates, primers = zip(*pairs)

        temps = melting_temps(templates, primers, 1e-6, 1e-6, 0, 100)

        sock.send(pickle.dumps(temps))


if __name__ == "__main__":

    import sys

    if sys.argv[1] == 'server':
        server()

    else:
        server_list = sys.argv[1:]
        randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))

        templates = [ randseq(40) for _ in range(24000) ]
        primers = [ randseq(40) for _ in range(24000) ]
        pairs = zip(templates, primers)

        print client(server_list, pairs)
