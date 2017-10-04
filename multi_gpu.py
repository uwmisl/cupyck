import zmq
import pickle
from ctypes import *
import numpy as np
import math

port = 2046

def client(server_list, sequences):

    nservers = len(server_list)
    nseqs = len(sequences)

    splits = range(0, nseqs + 1, nseqs / nservers)
    splits[-1] = nseqs
    ranges = zip(splits, splits[1:])

    chunks = [ sequences[start:end] for start, end in ranges ]

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
    lib = cdll.LoadLibrary("./pfunc-cuda.so")

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
        seqs = pickle.loads(msg)

        results = pfmulti(seqs, [37] * len(seqs), [1] * len(seqs))

        sock.send(pickle.dumps(results))


if __name__ == "__main__":

    import sys

    if sys.argv[1] == 'server':
        server()

    else:
        server_list = sys.argv[1:]
        randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))
        seqs = [ randseq(80) for _ in range(100000) ]
        print client(server_list, seqs)
