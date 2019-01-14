from ..cluster import Server
import numpy as np
import math

from nupyck.apps.concentrations import calc_conc

def pair_perms((template, primer)):
    return [ template,
             primer,
             "+".join((template, template)),
             "+".join((template, primer)),
             "+".join((primer, primer))
            ]

def melting_temps(session, templates, primers, t_conc, p_conc, t_lo, t_hi):
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
    Gs = np.array(
        session.pfunc(
            seqs.flatten(),
            temps,
             syms.flatten()
        )['energies']
    ).reshape(-1, 5)
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
    Gs = np.array(
        session.pfunc(
            seqs[rem].flatten(),
            temps,
            syms[rem].flatten()
        )['energies']
    ).reshape(-1, 5)
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
        Gs = np.array(
                session.pfunc(
                seqs[rem].flatten(),
                np.repeat(temps, 5),
                syms[rem].flatten()
            )['energies']
        ).reshape(-1,5)
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

class MtServer(Server):

    def __init__(self, t_conc, p_conc, t_lo, t_hi, **kwargs):

        kwargs['t_lo'] = t_lo
        kwargs['t_hi'] = t_hi

        super(MtServer, self).__init__(**kwargs)

        self.t_conc = t_conc
        self.p_conc = p_conc

        self.t_lo = t_lo
        self.t_hi = t_hi

    def worker(self, *args):
        templates, primers = args
        return melting_temps(
            self,
            templates,
            primers,
            self.t_conc, self.p_conc,
            self.t_lo, self.t_hi
        )


if __name__ == "__main__":

    import sys

    if sys.argv[1] == 'server':
        print "initializing..."
        server = MtServer(
            t_conc = 1e-8,
            p_conc = 1e-6,
            t_lo = 0,
            t_hi = 100,
            t_step = 0.05,
            port = 2046,
            nblocks = 16384,
            nthreads = 64,
            max_seqlen = 100
        )
        print "listening."
        server.listen()

    else:
        from ..cluster import Client
        import nupyck.apps.concentrations as nconc
        def check((template, primer)):
            return nconc.melting_temp(template, primer, 1e-8, 1e-6, 0, 100)

        server_list = sys.argv[1:]
        calc_melts = Client(server_list, 2046)

        randseq = lambda n: "".join(np.random.choice(list("ATCG"), n))
        templates = [ randseq(40) for _ in range(500) ]
        primers = [ randseq(40) for _ in range(500) ]
        pairs = zip(templates, primers)

        gpu_melts = np.array(calc_melts(pairs))
        cpu_melts = np.array(map(check, pairs))

        print np.mean(np.square(gpu_melts - cpu_melts))
