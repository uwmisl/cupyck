import zmq
import pickle
import abc
import sys
import time
import pandas as pd


class Server(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, port, session):

        ctxt = zmq.Context()
        self.sock = ctxt.socket(zmq.REP)
        self.sock.bind("tcp://*:%d" % port)

        self.session = session

    @abc.abstractmethod
    def worker(self, jobs):
        pass

    def listen(self, verbose = False):
        if verbose:
            print "listening..."

        while True:

            request = self.sock.recv()
            jobs = pickle.loads(request)

            if verbose:
                sys.stdout.write("processing %d jobs..." % len(jobs))
                start = time.time()

            result = self.worker(jobs)

            if verbose:
                end = time.time()
                sys.stdout.write("done in %f seconds.\n" % (end - start))

            reply = pickle.dumps(result)
            self.sock.send(reply)


class Client(object):

    def __init__(self, server_list):
        self.socks = []

        ctxt = zmq.Context()
        for server, port in server_list:
            sock = ctxt.socket(zmq.REQ)
            sock.connect("tcp://%s:%d" % (server, port))
            self.socks.append(sock)

    def __call__(self, items):

        nitems = len(items)
        nservers = len(self.socks)
        splits = range(0, nitems + 1, nitems / nservers)
        splits[-1] = nitems
        ranges = zip(splits, splits[1:])
        chunks = [ items[start:end] for start, end in ranges ]

        for idx, sock in enumerate(self.socks):
            request = pickle.dumps(chunks[idx])
            sock.send(request)

        results = []
        for sock in self.socks:
            reply = sock.recv()
            result = pickle.loads(reply)
            results.append(result)

        results = pd.concat(results)

        return results

    def close(self):
        for sock in self.socks:
            sock.close()
