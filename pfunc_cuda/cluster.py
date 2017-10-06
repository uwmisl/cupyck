from .session import Session
import zmq
import pickle
import itertools


class Server(Session):

    def __init__(self, port, *args, **kwargs):

        ctxt = zmq.Context()
        self.sock = ctxt.socket(zmq.REP)
        self.sock.bind("tcp://*:%d" % port)

        super(Server, self).__init__(*args, **kwargs)

    def worker(self, *args):
        return self.pfunc(*args)['energies']

    def listen(self):

        while True:

            request = self.sock.recv()
            items = pickle.loads(request)
            args = zip(*items)

            result = self.worker(*args)

            reply = pickle.dumps(result)
            self.sock.send(reply)


class Client:

    socks = []

    def __init__(self, server_list, port):

        ctxt = zmq.Context()
        for server in server_list:
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

        results = list(itertools.chain(*results))

        return results

