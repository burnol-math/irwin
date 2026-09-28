# -*- mode: python ; coding: utf-8; -*-

# irwin_v5_mpmath_spawn.py
# Python + mpmath conversion of irwin_v5.sage (no SageMath needed).

# THIS FILE WAS PRODUCED ON MONDAY SEP 28, 2026 BY CLAUDE OPUS 5.5
# AS A CONVERSION OF irwin_v5.sage, version 1.5.7 of 2025/05/17,
# The original was fully authored by Jean-François Burnol.
# THIS FILE IS A CONVERSION BY AN AI TO PYTHON + MPMATH + GMPY2
# PLUS USAGE OF SPAWN ON MACOS AS IT IS THE DEFAULT FOR PYTHON.

# THE HUMAN HAS NOT CHECKED ANYTHING, USE AT OWN RISK, INSTRUCTIONS
# BELOW. The same License apply as to the original, see below.

# For obvious deficiencies of the macOS Python check the author
# comments at gitlab.com/burnolmath/irwin.

# Use for example (interactive python or ipython):
#
#     >>> import irwin_v5_mpmath_spawn
#     >>> irwin_v5_mpmath_spawn.maxworkers = 10   # optional, default 8
#     >>> from irwin_v5_mpmath_spawn import irwin, irwinpos
#     >>> irwin(10, 9, 0, 52)
#     22.92067661926415034816365709437593191494476243699848
#
# In a script, the calls must be protected by
#
#     if __name__ == "__main__":
#
# because the worker processes are started with the "spawn" method
# (the macOS default), which re-imports the main module in each worker.

__version__  = "1.5.7-py"
__date__     = "2026/09/28"
__filename__ = "irwin_v5_mpmath_spawn.py"

irwin_v5_docstring = """
This file is a conversion to Python + mpmath of irwin_v5.sage,
version 1.5.7 of 2025/05/17, from

    https://gitlab.com/burnolmath/irwin

which is itself an evolution of irwin.sage as available at

    https://arxiv.org/src/2402.09083v5/anc

The SageMath RealField(p) objects are replaced by precisions p in
bits: every number is a raw mpmath "mpf" value (a tuple (sign,
mantissa, exponent, bitcount), see mpmath.libmp) and every
arithmetic operation is done with explicit rounding to nearest at
the precision which the SageMath code would have used.  Having
gmpy2 installed is strongly recommended: mpmath then uses GMP for
its integer arithmetic.

The SageMath @parallel decorator (which forks at each call) is
replaced by maxworkers persistent worker processes, started with
the "spawn" method, which keep their own copies of the data they
need and receive only the new coefficients u_{j;m} or v_{j;m} at
each call.  Each worker also maintains its own row of the Pascal
triangle.

Main features of irwin_v5.sage (kept here):

- Use of decreasing precision for higher terms in the series.
  This is done with a granularity of PrecStep bits.  It is an
  optional parameter for the procedures irwin() and irwinpos().
  It defaults to 500.

- There is no pre-computation of the first 1000 rows of the Pascal
  triangle.  Only maxworkers rows of the Pascal triangle are kept
  at any given time in memory in the main process, and each worker
  keeps one row.  See file taille_pascal.pdf for details on the
  storage size needed for rows of the Pascal triangle.

- Parallelization with maxworkers processes, where maxworkers
  defaults to 8 and can be modified via irwin_v5_mpmath_spawn.maxworkers = N
  (no attempt is made to adjust dynamically to the number of cores
  available).  The new value is taken into account at the next call
  of irwin() or irwinpos(), there is no need to reload the module.

Miscellaneous remarks (from irwin_v5.sage):

* The default level is 3.

* All printed values round away the nbguardbits guard bits used
  internally for the computations.

* The global variable nbguardbits has initial value 12.

Original SageMath code:
Copyright (C) 2025, 2026 Jean-François Burnol
License: CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/
This Python conversion is distributed under the same license.
"""

import atexit
import math
import multiprocessing as _mp
import pickle
import time
import traceback

import mpmath
from mpmath.libmp import (MPZ, from_int, fzero, fone,
                          mpf_add, mpf_mul, mpf_div, mpf_neg, mpf_pos,
                          to_str, to_float)

_RND = 'n'  # round to nearest, as SageMath RealField default


def _fillin_irwin_docstring():
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(
            irwin_v5_fn_docstring.format('u_{j;m}', 'irwin') if
            (obj.__name__ == 'irwin') else
            irwin_v5_fn_docstring.format('v_{j;m}', 'irwinpos')
        )
        return obj
    return dec

irwin_v5_fn_docstring = """:param int b: the integer base
    :param int d: the digit.
    :param int k: the number of occurrences.
    :param int nbdigits: (optional, default 34)
        The wished-for number of decimal digits for the result.
    :param int level: (optional, default 3)
        The level must be 2, 3 or 4.
    :param int PrecStep: (optional, default ``500``)
        Terms of the series are computed with an evolving
        precision, which differs from the maximal precision by a
        suitable multiple of PrecStep.
    :param bool all: (optional, default ``False``)
        If ``True``  print all irwin sums for ``j`` occurrences
        with ``j`` from ``0`` to ``k``.
    :param bool showtimes: (optional, default ``False``)
        Whether to print out timings for various steps.
    :param bool verbose: (optional, default ``False``)
        Whether to print the values of intermediate contributions
        to the final value, in particular to confirm enough terms
        of the series were used.
    :param int Mmax: (optional, default ``-1``)
        If not ``-1`` the number of terms of the series to use.
        Use only if the auto-choice is excessive (this will be the
        case for k=0 and d=1, and to a lesser extent when d=b-1).
        Use ``verbose=True`` to check how many terms are used by
        default.
    :param bool persistentpara: (optional, default ``True``)
        Whether once parallel mode is chosen to compute the
        coefficients {0}'s to check again if non-parallel would be
        better.

    :rtype: :class:`IrwinNumber` (a subclass of mpmath.mpf)
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b.

    The returned value is an mpmath mpf holding the result rounded
    to the final binary precision; it prints with nbdigits
    significant decimal digits.

    The Burnol algorithm depends on a choice of "level".  The
    default is level=3 which is appropriate for obtaining
    hundreds of digits or more.  Setting level=4 seems to be
    useful only for small bases b such as b=2 or 3, it seems not
    to be useful for b=10.  The higher the k parameter is
    (required number of occurrences of the digit d), the sooner
    level=3 (which is default) is better choice than level=2.
    Actual thresholds may depend on the maxworkers setting and
    number of cores actually available on your computing device,
    which the code does not try to query, it is up to user to
    set irwin_v5_mpmath_spawn.maxworkers appropriately.

    Example:
    --------

    >>> {1}(10, 9, 4, 52, all=True)
    (k=0) 22.92067661926415034816365709437593191494476243699848
    (k=1) 23.04428708074784831967594930973617482538959203064774
    (k=2) 23.02604026596124378845022249787272342108112267542086
    (k=3) 23.02585299837244431714290384468012275518705238435290
    (k=4) 23.02585095265829261377053973815542996035002267989413
    23.02585095265829261377053973815542996035002267989413
"""

nbguardbits = 12
maxworkers = 8
maxworkersinfostring = ("maxworkers variable has been created "
                        "and assigned value 8.")


# ----------------------------------------------------------------------
# Exact integer helpers replacing SageMath exact symbolic floor/ceil.
# ----------------------------------------------------------------------

def _max_pow_le(b, e):
    """Largest integer t >= 0 with b**t <= 2**e (b >= 2, e >= 0)."""
    B = 1 << e
    t = int(e / math.log2(b))
    while t > 0 and b ** t > B:
        t -= 1
    while b ** (t + 1) <= B:
        t += 1
    return t


def _min_pow_ge(b, e):
    """Smallest integer s >= 0 with b**s >= 2**e (b >= 2, e >= 0)."""
    t = _max_pow_le(b, e)
    return t if b ** t == (1 << e) else t + 1


def _nbdigits_for_prec(prec):
    """Number of decimal digits used to print a value of precision prec.

    It is floor(prec*log10(2)) - 1 (at least 1).  With the choice of
    nbbits_final = ceil((nbdigits+1)*log2(10)) this gives back
    nbdigits, as in the SageMath outputs.
    """
    return max(1, _max_pow_le(10, prec) - 1)


# ----------------------------------------------------------------------
# Emulation of RealField(p) coercions and printing.
# ----------------------------------------------------------------------

def _R(x, p):
    """Emulates RealField(p)(x) for x a Python int or an mpf tuple."""
    if isinstance(x, tuple):
        return mpf_pos(x, p, _RND)
    return from_int(MPZ(x), p, _RND)


def _plus(x, y, p):
    """x + y where x, y are ints (Sage Integers) or mpf tuples."""
    if not isinstance(x, tuple) and not isinstance(y, tuple):
        return x + y
    return mpf_add(_R(x, p), _R(y, p), p, _RND)


def _times(x, y, p):
    """x * y where x, y are ints (Sage Integers) or mpf tuples."""
    if not isinstance(x, tuple) and not isinstance(y, tuple):
        return x * y
    return mpf_mul(_R(x, p), _R(y, p), p, _RND)


def _suminv(L, p):
    """sum(1/RealField(p)(x) for x in L), which is int 0 if L is empty."""
    S = 0
    for x in L:
        S = _plus(S, mpf_div(fone, from_int(MPZ(x), p, _RND), p, _RND), p)
    return S


def _fmt(x, prec):
    """String for printing x as SageMath would print a RealNumber."""
    if not isinstance(x, tuple):
        return str(x)
    return to_str(x, _nbdigits_for_prec(prec), strip_zeros=False)


class IrwinNumber(mpmath.mpf):
    """An mpmath mpf which prints with a given number of digits."""
    _irwin_nbdigits = 15

    def __str__(self):
        return to_str(self._mpf_, self._irwin_nbdigits, strip_zeros=False)

    __repr__ = __str__

    def __format__(self, spec):
        # f-strings call __format__, which recent mpmath versions define
        # for mpf and which would otherwise ignore __str__ above.
        if not spec:
            return str(self)
        return mpmath.mpf.__format__(self, spec)


def _make_result(x, prec):
    """Emulates returning Rfinal(S)."""
    r = IrwinNumber()
    # set the value directly: IrwinNumber(x) would round x to the
    # current mpmath working precision mpmath.mp.prec
    r._mpf_ = mpf_pos(x, prec, _RND)
    r._irwin_nbdigits = _nbdigits_for_prec(prec)
    return r


# ----------------------------------------------------------------------
# Partial recurrences for the u_{j;m} or v_{j;m} (shared by main
# process in serial mode and by the workers in parallel mode).
# ----------------------------------------------------------------------

def _v5_ukm_partial(a, m, Pm, G, D, T, p, k):
    """Recurrences (partial) for the u_{j;m}'s or v_{j;m}'s.

    - This handles all j's from 0 to k (because to compute
      for k we need for k-1, and to compute for k-1, we
      need for k-2 and so on until j=0).
    - Pm[i] stands for binomial coefficient "m choose i".
    - G[i] stands for gamma (or gammaprime) power sum.
    - D[i] is d**i or dprime**i (dprime = b-1-d)
    - T[n][j] holds previously known u_{j;n}'s or v_{j;n}'s.
    - p is the precision (in bits) suitable for the evaluation
      of the u_{j;m}'s or v_{j;m}'s (the RealField Rm of the
      SageMath code).

    The formulas are those of Theorem 1 (equations (2) and (3))
    and Theorem 4 (equations (5) and (6)) in the numeration as
    in arXiv:2402.09083v5.  Up to the division by b**(m+1)-b+1
    which will be done later, because finitely many terms are
    still lacking at this stage.

    When the procedure is called we have computed all u_{j;n}'s
    or v_{j;n}'s up to some M.  The procedure will be called in
    parallel with m=M+1, M+2, ...., M+n with a=1, 2, ..., n.
    The quantity m-a is thus the same M for this parallelized
    bunch.  Once this procedure returns we have value for m=M+1
    exactly, but will need correction for M+2, then M+3, ... up
    to the last one M+n where n is most of the time a multiple
    of maxworkers.  And there will be division by b**(m+1)-b+1.

    Memo: for j=0 and the v_{0;m}'s there is an extra contribution
    b**(m+1) which is added by the caller.
    """
    # Pm[i]*Rm(G[i]) does not depend on j, compute it only once.
    PG = [mpf_mul(from_int(Pm[i], p, _RND), mpf_pos(G[i], p, _RND),
                  p, _RND) for i in range(a, m + 1)]
    if k > 0:
        PD = [mpf_mul(from_int(Pm[i], p, _RND), mpf_pos(D[i], p, _RND),
                      p, _RND) for i in range(a, m + 1)]
    result = []
    for j in range(k + 1):
        A = fzero
        for i in range(a, m + 1):
            A = mpf_add(A, mpf_mul(PG[i - a], mpf_pos(T[m - i][j], p, _RND),
                                   p, _RND), p, _RND)
        if j >= 1:
            B = fzero
            for i in range(a, m + 1):
                B = mpf_add(B, mpf_mul(PD[i - a],
                                       mpf_pos(T[m - i][j - 1], p, _RND),
                                       p, _RND), p, _RND)
            A = mpf_add(A, B, p, _RND)
        result.append(A)
    return result


# ----------------------------------------------------------------------
# Worker processes.
# ----------------------------------------------------------------------

def _pascal_row(m):
    """Row m of the Pascal triangle, computed from scratch."""
    row = [MPZ(1)] * (m + 1)
    c = MPZ(1)
    for i in range(1, m // 2 + 1):
        c = c * (m - i + 1) // i
        row[i] = c
        row[m - i] = c
    return row


def _worker_setup(st, precs, A1, dval, k):
    st.clear()
    st['precs'] = precs
    st['k'] = k
    st['cur'] = precs[0]
    st['T'] = []
    st['A1'] = [MPZ(a) for a in A1]
    st['pw'] = [MPZ(1) for a in A1]    # a**e for a in A1
    st['G'] = [fzero]                  # G[0] is never used
    st['dval'] = MPZ(dval)
    st['dpw'] = MPZ(1)                 # dval**e
    st['D'] = [fone]                   # D[0] is never used
    st['prow'] = None
    st['pm'] = -1


def _worker_ukm(st, a, M, start, rows):
    """Compute _v5_ukm_partial for m = M + a.

    rows is the pickled list of the coefficients rows start, ..., M
    which the worker does not already have.  The worker keeps them rounded to the
    precision needed for m = M+1 (this is at least the precision
    needed for all later m's), and rounds again when the precision
    decreases, to keep its memory footprint as low as possible.
    """
    precs = st['precs']
    k = st['k']
    T = st['T']
    G = st['G']
    D = st['D']
    if len(T) != start:
        raise RuntimeError(f"worker has {len(T)} rows, expected {start}")
    rows = pickle.loads(rows)
    m = M + a
    tp = precs[M + 1]
    if tp < st['cur']:
        st['cur'] = tp
        for n in range(len(T)):
            T[n] = [mpf_pos(x, tp, _RND) for x in T[n]]
        for i in range(1, len(G)):
            G[i] = mpf_pos(G[i], tp, _RND)
        if k > 0:
            for i in range(1, len(D)):
                D[i] = mpf_pos(D[i], tp, _RND)
    cur = st['cur']
    for row in rows:
        T.append([mpf_pos(x, cur, _RND) for x in row])

    # gammas: G[i] = RealField(precs[i])(sum(a**i for a in A1))
    A1 = st['A1']
    pw = st['pw']
    while len(G) <= m:
        i = len(G)
        for n in range(len(A1)):
            pw[n] *= A1[n]
        g = from_int(sum(pw, MPZ(0)), precs[i], _RND)
        G.append(mpf_pos(g, cur, _RND))
    if k > 0:
        while len(D) <= m:
            i = len(D)
            st['dpw'] *= st['dval']
            D.append(mpf_pos(from_int(st['dpw'], precs[i], _RND),
                             cur, _RND))

    # Pascal triangle row m
    pm = st['pm']
    if pm == m:
        prow = st['prow']
    elif 0 <= pm < m and m - pm <= 64:
        prow = st['prow']
        for _ in range(m - pm):
            n = len(prow)
            prow = ([prow[0]] + [prow[i - 1] + prow[i] for i in range(1, n)]
                    + [prow[0]])
    else:
        prow = _pascal_row(m)
    st['prow'] = prow
    st['pm'] = m

    return _v5_ukm_partial(a, m, prow, G, D, T, precs[m], k)


def _worker_beta(st, start, end, step, nblock):
    """For m in range(start, end, step) sum of 1/n**(m+1), n in nblock.

    Each 1/n**(m+1) is 1/RealField(precs[m])(n**(m+1)), as in the
    SageMath code.
    """
    precs = st['precs']
    res = []
    if start >= end:
        return res
    powers = [MPZ(n) ** (start + 1) for n in nblock]
    nstep = [MPZ(n) ** step for n in nblock]
    m = start
    while True:
        p = precs[m]
        S = fzero
        for x in powers:
            S = mpf_add(S, mpf_div(fone, from_int(x, p, _RND), p, _RND),
                        p, _RND)
        res.append(S)
        m += step
        if m >= end:
            break
        powers = [x * y for x, y in zip(powers, nstep)]
    return res


def _worker_main(conn):
    """Loop of a worker process."""
    st = {}
    while True:
        try:
            msg = conn.recv()
        except EOFError:
            break
        op = msg[0]
        if op == "quit":
            break
        try:
            if op == "setup":
                _worker_setup(st, *msg[1:])
                res = None
            elif op == "clear":
                st.clear()
                res = None
            elif op == "ukm":
                res = _worker_ukm(st, *msg[1:])
            elif op == "beta":
                res = _worker_beta(st, *msg[1:])
            else:
                raise ValueError(f"unknown operation {op}")
            conn.send((True, res))
        except Exception:
            conn.send((False, traceback.format_exc()))
    conn.close()


class _Workers:
    """maxworkers persistent processes, one Pipe each."""

    def __init__(self, n):
        ctx = _mp.get_context("spawn")
        self.n = n
        self.conns = []
        self.procs = []
        for _ in range(n):
            parent, child = ctx.Pipe()
            proc = ctx.Process(target=_worker_main, args=(child,),
                               daemon=True)
            proc.start()
            child.close()
            self.conns.append(parent)
            self.procs.append(proc)

    def alive(self):
        return all(proc.is_alive() for proc in self.procs)

    @staticmethod
    def _check(reply):
        ok, res = reply
        if not ok:
            raise RuntimeError(f"worker process raised:\n{res}")
        return res

    def send_bytes(self, w, data):
        self.conns[w].send_bytes(data)

    def recv(self, w):
        return self._check(self.conns[w].recv())

    def broadcast(self, msg):
        data = pickle.dumps(msg, protocol=pickle.HIGHEST_PROTOCOL)
        for w in range(self.n):
            self.conns[w].send_bytes(data)
        for w in range(self.n):
            self.recv(w)

    def close(self):
        for conn in self.conns:
            try:
                conn.send(("quit",))
                conn.close()
            except Exception:
                pass
        for proc in self.procs:
            proc.join(timeout=5)
            if proc.is_alive():
                proc.terminate()
        self.conns = []
        self.procs = []


_WORKERS = None


def _get_workers():
    """Return the pool of workers, (re)starting it if needed."""
    global _WORKERS
    if (_WORKERS is None or _WORKERS.n != maxworkers
            or not _WORKERS.alive()):
        if _WORKERS is not None:
            _WORKERS.close()
        _WORKERS = _Workers(maxworkers)
    return _WORKERS


def shutdown_workers():
    """Terminate the worker processes (they are restarted when needed)."""
    global _WORKERS
    if _WORKERS is not None:
        _WORKERS.close()
        _WORKERS = None


atexit.register(shutdown_workers)


# ----------------------------------------------------------------------
# Recurrences for the u_{j;m}'s or v_{j;m}'s.
# ----------------------------------------------------------------------

def _v5_umtimeinfo(single_ns, multi_ns, para, wrkrs, M, s):
    """Auxiliary shared between irwin() and irwinpos().
    """
    multi = multi_ns * 1e-9
    single = single_ns * 1e-9
    if multi_ns < single_ns:
        if para:
            print(f"... poursuite car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle ({M}<m<={M+s})")
        else:
            print(f"... basculement car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle (maxworkers={wrkrs}, {M}<m<={M+s})")
    else:
        if para:
            print(f"... on quitte car {single:.3f}s<{multi:.3f}s"
                  f" l'exécution parallèle ({M}<m<={M+s})")
        else:
            print(f"... pas utile ({single:.3f}s<{multi:.3f}s)"
                  f" d'exécuter en parallèle ({M}<m<={M+s})")


def _v5_setup_para_recurrence(touslescoeffs, Gammas, PuissancesDeD,
                              PascalRows, IndexToPrec, b, bmoinsun, k,
                              showtimes, persistentpara, is_for_vm,
                              workers):
    """Set up procedure calling _v5_ukm_partial and completing its job.
    """
    nworkers = workers.n
    # synced[w] = number of rows of touslescoeffs known to worker w
    synced = [0] * nworkers

    def _parallel_partial(M, step):
        """Let worker a-1 compute _v5_ukm_partial for m=M+a, a=1..step.

        Each worker first receives the rows of touslescoeffs it does
        not have yet.  The pickling of these rows is done only once
        for all workers which need the same rows.
        """
        cache = {}
        for a in range(1, step + 1):
            w = a - 1
            start = synced[w]
            if start not in cache:
                cache[start] = pickle.dumps(touslescoeffs[start:M + 1],
                                            protocol=pickle.HIGHEST_PROTOCOL)
            data = pickle.dumps(("ukm", a, M, start, cache[start]),
                                protocol=pickle.HIGHEST_PROTOCOL)
            workers.send_bytes(w, data)
            synced[w] = M + 1
        ukm_partial = [None]
        for a in range(1, step + 1):
            ukm_partial.append(workers.recv(a - 1))
        return ukm_partial

    def _serial_partial(M, step):
        m = M
        ukm_partial = [None]
        for j in range(1, 1 + step):
            m += 1
            ukm_partial.append(_v5_ukm_partial(j, m, PascalRows[j],
                                               Gammas, PuissancesDeD,
                                               touslescoeffs,
                                               IndexToPrec[m], k))
        return ukm_partial

    def _v5_para_recurrence(m, step, useparallel):
        """Wrapper of parallelized calls to _v5_ukm_partial().

        First we compute "step" (which is maxworkers or less than it)
        new rows of the Pascal triangle of binomial coefficients.  We
        use some specificities of how Python handles list type to do
        that in a way persistent in memory across calls.

        Then, if useparallel is True we let the workers compute
        _v5_ukm_partial() for m varying from M+1 to M+step, where M
        is the initial value of argument m.  If useparallel is False
        we still do that from time to time to compare with computing
        serially.  Even with useparallel True and except if
        persistentpara is False we will check from time to time the
        comparison between parallel and serial.

        If useparallel is False, we compute serially new u_{j;m}'s
        or v_{j;m}'s.

        In all cases the formulas of arXiv:2402.09083 are applied.
        In order to share code, when computing serially we do as in
        the parallel branch and first evaluate only partially the
        recurrent formulas.  So we can then correct via the missing
        terms in both cases.  The recurrences for the v_{j;m}'s
        differ from those for the u_{j;m}'s in what Gammas and
        PuissancesDeD stand for, as well as one unique extra term in
        the recurrence computing v_{0;m}.

        We do the computation of the u_{j;m}, v_{j;m} for given m
        from j=0 upto j=k.  This gives a list which is appended to
        the list touslescoeffs holding all such coefficients (so the
        indexing is in reverse order compared to the mathematical
        notation: m first, and j second).
        """
        del PascalRows[:-1]
        for i in range(step):
            m += 1
            newPascalRow = [ 1 ]
            halfm = m // 2
            newPascalRow.extend([PascalRows[-1][j-1] +
                                 PascalRows[-1][j] for j in range(1, halfm)])
            if not (m&1):
                halfPascalRow = newPascalRow.copy()
            newPascalRow.append(PascalRows[-1][halfm-1]+PascalRows[-1][halfm])
            if m&1:
                newPascalRow.extend(reversed(newPascalRow))
            else:
                newPascalRow.extend(reversed(halfPascalRow))
            if i == 0:
                del PascalRows[:]
                PascalRows.append(None)
            PascalRows.append(newPascalRow)

        M = m - step
        if ((M - 400) % 500 < nworkers):
            starttime_ns = time.perf_counter_ns()
            ukm_partial = _parallel_partial(M, step)
            multitime_ns = time.perf_counter_ns() - starttime_ns

            if useparallel and persistentpara:
                # do not check again
                if showtimes:
                    print(f"... mode parallèle persistant "
                          f"({multitime_ns*1e-9:.3f}s; "
                          f"{M}<m<={M+step})")
            else:
                starttime_ns = time.perf_counter_ns()
                ukm_partial = _serial_partial(M, step)
                singletime_ns = time.perf_counter_ns() - starttime_ns

                if showtimes:
                    _v5_umtimeinfo(singletime_ns, multitime_ns,
                                   useparallel,
                                   nworkers, M, step)

                useparallel = multitime_ns < singletime_ns

        elif useparallel:
            ukm_partial = _parallel_partial(M, step)

        else:
            ukm_partial = _serial_partial(M, step)

        # Now correct the um's (or vm's) (prior to dividing by b**(m+1)-b+1)
        # via the addition of finitely missing contributions in order of
        # increasing m's.
        m = M
        for j in range(1, 1 + step):
            m += 1
            p = IndexToPrec[m]
            Pj = PascalRows[j]
            D = from_int(MPZ(b**(m+1) - bmoinsun), p, _RND)
            # Attention to the b**(m+1) extra term specific to v_m recurrence.
            if is_for_vm:
                num = mpf_add(from_int(MPZ(b ** (m+1)), p, _RND),
                              ukm_partial[j][0], p, _RND)
            else:
                num = mpf_pos(ukm_partial[j][0], p, _RND)
            s = fzero
            for i in range(1, j):
                s = mpf_add(s, mpf_mul(mpf_mul(from_int(Pj[i], p, _RND),
                                               mpf_pos(Gammas[i], p, _RND),
                                               p, _RND),
                                       mpf_pos(touslescoeffs[m-i][0],
                                               p, _RND),
                                       p, _RND), p, _RND)
            cm = [ mpf_div(mpf_add(num, s, p, _RND), D, p, _RND) ]
            for q in range(1, k+1):
                s1 = fzero
                for i in range(1, j):
                    s1 = mpf_add(s1,
                                 mpf_mul(mpf_mul(from_int(Pj[i], p, _RND),
                                                 mpf_pos(Gammas[i], p, _RND),
                                                 p, _RND),
                                         mpf_pos(touslescoeffs[m-i][q],
                                                 p, _RND),
                                         p, _RND), p, _RND)
                s2 = fzero
                for i in range(1, j):
                    s2 = mpf_add(s2,
                                 mpf_mul(mpf_mul(from_int(Pj[i], p, _RND),
                                                 mpf_pos(PuissancesDeD[i],
                                                         p, _RND),
                                                 p, _RND),
                                         mpf_pos(touslescoeffs[m-i][q-1],
                                                 p, _RND),
                                         p, _RND), p, _RND)
                num = mpf_add(ukm_partial[j][q], s1, p, _RND)
                num = mpf_add(num, cm[-1], p, _RND)
                num = mpf_add(num, s2, p, _RND)
                cm.append(mpf_div(num, D, p, _RND))
            touslescoeffs.append(cm)
        # Update status.
        return useparallel
    return _v5_para_recurrence


# ----------------------------------------------------------------------
# The beta's (sums of inverse powers).
# ----------------------------------------------------------------------

def _v5_beta(workers, inputdata):
    """Parallelized computation of beta coefficients.

    inputdata is a list of (start, end, nblock) with at most
    maxworkers entries.  For each entry, m varies in
    range(start, end, maxworkers) and we compute the sum of
    1/n**(m+1) for n varying in given "nblock".  IndexToPrec maps
    indices m to suitable precision.  Higher indices use lower
    precision, this is why the range of m's is split according to
    value modulo maxworkers, so that the computation costs are
    about equal across workers.

    Returns the list of the results, in the order of inputdata.
    """
    for w, (start, end, nblock) in enumerate(inputdata):
        workers.send_bytes(w, pickle.dumps(("beta", start, end,
                                            workers.n, nblock),
                                           protocol=pickle.HIGHEST_PROTOCOL))
    return [workers.recv(w) for w in range(len(inputdata))]


def _v5_map_beta_notimes(Mmax, workers, maxblock):
    """Sets up a procedure to call _v5_beta() and assembles its results.

    The defined procedure will receive an argument j which is in the
    range from 0 to k inclusive.  It will then use the integers in
    the "block" maxblock[j] as the ones for which the sum of inverse
    powers needs to be computed.

    The procedure defined by this does not display intermediate
    computing times.
    """
    W = workers.n
    def map__v5_beta(j):
        """Calls parallelized _v5_beta() and assembles its results.

        After having computed beta_{m+1}'s for m's split by
        their modulo maxworkers value (in (1,..., maxworkers))
        we reorganize the maxworkers lists of values into a
        single list L in order of increasing m's, L[0] = 0.
        """
        L = [0] * (Mmax + 1)
        inputdata = [(i, Mmax + 1, maxblock[j]) for i in range(1, 1 + W)]
        results_1 = _v5_beta(workers, inputdata)
        for i, res in zip(range(1, 1 + W), results_1):
            L[i:Mmax + 1:W] = res
        return L
    return map__v5_beta


def _v5_map_beta_withtimes(Mmax, workers, maxblock):
    """Sets up a procedure to call _v5_beta() and assembles its results.

    The defined procedure will receive an argument j which is in the
    range from 0 to k inclusive.  It will then use the integers in
    the "block" maxblock[j] as the ones for which the sum of inverse
    powers needs to be computed.

    The procedure defined by this does displays intermediate
    computing times.  It will divide the range from 1 to
    Mmax in chunks of size a multiple of maxworkers near to 1000.
    If maxworkers if 32 or more, chunks of size 32*maxworkers are
    used for displaying their timings.
    """
    W = workers.n
    # We want to display some visual sign of progress.
    # Find the largest multiple of maxworkers at most 1000,
    # do something reasonable if maxworkers is big
    q = max(1000 // W, 32)
    mSize = q * W
    def map__v5_beta(j):
        """Calls parallelized _v5_beta() and assembles its results.

        And compute intermediate timings while doing it.
        """
        print(f"... ({j} occ.) ", end = "", flush = True)
        starttime = time.perf_counter()
        lasttime = starttime
        mbegin = 1  # will remain congruent to 1 modulo maxworkers
        mend = 1    # this one also
        L = [0] * (Mmax + 1)
        for rep in range(Mmax // mSize):
            mend   = mbegin + mSize
            inputdata = [(mbegin + i, mend, maxblock[j]) for i in range(W)]
            results_1 = _v5_beta(workers, inputdata)
            for i, res in zip(range(W), results_1):
                L[mbegin + i:mend:W] = res
            stoptime = time.perf_counter()
            print(f"m<{mend} ({stoptime-lasttime:.3f}s)",
                  end = "\n             ", flush= True)
            lasttime = stoptime
            mbegin = mend

        if mend < Mmax+1:
            inputdata = [(mend + i, Mmax + 1, maxblock[j]) for i in range(W)]
            results_1 = _v5_beta(workers, inputdata)
            for i, res in zip(range(W), results_1):
                L[mend + i:Mmax + 1:W] = res
        stoptime = time.perf_counter()
        if mend < Mmax + 1:
            print(f"m<{Mmax+1} ({stoptime-lasttime:.3f}s)",
                  end = " ")
        print(f"Fini! En tout : {stoptime-starttime:.3f}s")
        return L
    return map__v5_beta


def _v5_shorten_small_real(rr):
    """Get magnitude order of a tiny real number

    Probably very clumsy.
    """
    with mpmath.workprec(53):
        x = mpmath.mpf(mpf_pos(rr, 53, _RND))
        s = 1 if x >= 0 else -1
        y = mpmath.log10(abs(x))
        N = int(y)  # truncation towards zero
        return s*10**float(y-N+1), N-1


def _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax=-1):
    """Preparation of an array mapping each m to a precision.

    See irwin_v5_doc.pdf for mathematical details.

    All floor() and ceil() of the SageMath code are computed
    exactly via integer arithmetic.
    """

    # Chose number of bits to (try to) guarantee we will have nbdigits
    # decimal digits in output.
    # nbbits_final = ceil((nbdigits+1)*log(10,2))
    nbbits_final = (10 ** (nbdigits + 1)).bit_length()

    # Computations are done (for the main terms) with elevated precision.
    nbbits = nbbits_final + nbguardbits

    # See irwin_v5_doc.pdf for the mathematical justification for this
    # choice of Mmax, which is the number of terms used from the
    # series given in Burnol papers.
    # _Mmax = floor((nbbits - nbguardbits/2)/(level-1)/log(b,2))
    # is the largest M with b**(2*M*(level-1)) <= 2**(2*nbbits-nbguardbits).
    _Mmax = _max_pow_le(b, 2 * nbbits - nbguardbits) // (2 * (level - 1))

    # The number of distinct precisions we need.
    # NbOfPrec = 1 + floor((nbbits - nbguardbits/2)/PrecStep)
    NbOfPrec = 1 + (2 * nbbits - nbguardbits) // (2 * PrecStep)
    LesPrecs = [nbbits - j * PrecStep for j in range(NbOfPrec)]

    # See irwin_v5_doc.pdf for the justification that we only need
    #
    #    nbbits - (l-1) * m * log(b,2)
    #
    # precision for the computation of the mth term.  So j is
    # chosen to be the largest such that nbbits - jT is at least
    # that value.  Hence j is floor((l-1)*log(b,2)*m/T) (with T =
    # PrecStep).
    #
    # We need to know which m will use a given j.
    # They verify j<= (l-1)*log(b,2)*m/T < j+1, i.e.
    #
    # ceil(j*T/(l-1)/log(b,2))<= m < ceil((j+1)*T/(l-1)/log(b,2)).
    #
    # ceil((j+1)*T/(l-1)/log(b,2)) is the smallest M such that
    # b**(M*(l-1)) >= 2**((j+1)*T), it is ceil(s/(l-1)) with s the
    # smallest integer such that b**s >= 2**((j+1)*T).
    IndexToPrec = []
    oldindexbound = 0
    for j in range(NbOfPrec):
        s = _min_pow_ge(b, (j + 1) * PrecStep)
        newindexbound = -((-s) // (level - 1))
        IndexToPrec.extend([LesPrecs[j]] * (newindexbound - oldindexbound))
        oldindexbound = newindexbound

    # We want m=1 to use always maximal precision.
    IndexToPrec[1] = nbbits
    IndexToPrec.extend([LesPrecs[-1]] * (_Mmax + 1 - len(IndexToPrec)))

    if Mmax == -1:
        Mmax = _Mmax
    else:
        # If user specified a custom Mmax which is beyond our
        # estimate we need to extend.
        if Mmax > _Mmax:
            IndexToPrec.extend([LesPrecs[-1]] * (Mmax - _Mmax))
            print(f"!!!! Warning Mmax={Mmax} is probably needlessly big.")
            print(f"!!!! {_Mmax} should be enough but we will use your value.")
            print(f"!!!! Use verbose=True to see the size of the smallest term.")

    return nbbits, nbbits_final, Mmax, IndexToPrec, NbOfPrec


def _v5_setup_blocks(b, d, level):
    """Organize integers according to nb of digits and d-count.

    It returns a list of lists of lists: blocks[l][j] is the list
    of integers having (l+1) digits in radix b, among whose exactly j
    are equal to d.  The "level" parameter is the maximal "l+1".

    - In particular for l=0, length-1 integers are exactly the non
    zero digits. blocks[0] always contains 2 entries
      * first one is the list of all non-zero digits distinct from d,
      * second one is either [d] or [] whether d is non-zero or zero.
    - blocks[1][0] = list of 2-digits integers (NOT strings!) with no d.
      blocks[1][1] = list of 2-digits integers with one occurrence of d.
      blocks[1][2] = [b*d+d] if d is not zero else [].
    - blocks[2][0] = list of 3-digits integers all whose digits are distinct
                     from d.
      blocks[2][1] = list of 3-digits integers with one digit equal to d.
      blocks[2][2] = list of 3-digits integers with two digits equal to d.
      blocks[2][3] = [b*b*d+b*d+d] if d is not zero else [].
    - idem for blocks[3] regarding 4-digits integers.
    We stop there as level accepted values are only 2, 3 or 4.
    """
    # A is the list of digits (inclusive of 0) not equal to d.
    A = [i for i in range(b)]
    A.remove(d)

    blocks = []

    block1 = []
    # block1[0]: pas le chiffre d (était aussi noté A1 dans code 2024)
    block1.append([a for a in A if a !=0])
    if d == 0:
        block1.append([])
    else:
        block1.append([d])

    blocks.append(block1)

    block2 = []
    # block2[0] = nombres à deux chiffres sans le chiffre d
    # block2[1] = nombres à deux chiffres avec une occurrence de d
    # block2[2] = nombres à deux chiffres avec deux occurrences de d
    #             sera présent mais vide si d=0
    block2.append([b * x + a for x in block1[0] for a in A])  # k=0
    L = [b * x + d for x in block1[0]]
    L.extend([b * x + a for x in block1[1] for a in A])
    block2.append(L)  # k = 1
    if d == 0:
        block2.append([])  # k = 2
    else:
        block2.append([(b + 1) * d])  # k = 2

    blocks.append(block2)

    if level > 2:
        block3 = []
        block3.append([b * x + a for x in block2[0] for a in A])  # k=0
        L = [b * x + d for x in block2[0]]
        L.extend([b * x + a for x in block2[1] for a in A])
        block3.append(L)  # k=1
        L = [b * x + d for x in block2[1]]
        L.extend([b * x + a for x in block2[2] for a in A])
        block3.append(L)  # k=2
        if d == 0:
            block3.append([])  # no length 3 number ddd if d=0
        else:
            block3.append([(b*b + b + 1) * d])  # k = 3

        blocks.append(block3)

    if level > 3:
        block4 = []
        block4.append([b * x + a for x in block3[0] for a in A])  # k=0
        L = [b * x + d for x in block3[0]]
        L.extend([b * x + a for x in block3[1] for a in A])
        block4.append(L)  # k=1
        L = [b * x + d for x in block3[1]]
        L.extend([b * x + a for x in block3[2] for a in A])
        block4.append(L)  # k=2
        L = [b * x + d for x in block3[2]]
        L.extend([b * x + a for x in block3[3] for a in A])
        block4.append(L)  # k=3

        if d == 0:
            block4.append([])
        else:
            block4.append([(b*b*b + b*b + b + 1) * d])  # k = 4

        blocks.append(block4)

    return blocks


def _v5_gammas(A1, bmoinsun, Mmax, IndexToPrec):
    """lesgammas[m] = RealField(IndexToPrec[m])(sum(a**m for a in A1))."""
    L = [ bmoinsun ]
    A = [MPZ(a) for a in A1]
    pw = [MPZ(1) for a in A1]
    for m in range(1, Mmax+1):
        for n in range(len(A)):
            pw[n] *= A[n]
        L.append(from_int(sum(pw, MPZ(0)), IndexToPrec[m], _RND))
    return L


def _v5_powers(dval, Mmax, IndexToPrec):
    """[1, RealField(IndexToPrec[m])(dval**m) for m in 1..Mmax]."""
    L = [ 1 ]
    x = MPZ(1)
    for m in range(1, Mmax+1):
        x *= dval
        L.append(from_int(x, IndexToPrec[m], _RND))
    return L


def _v5_series_term(touslescoeffs, betas, m, j, level, p):
    """Sum over i of touslescoeffs[m][j-i] * betas[i][m].

    betas[i] is the list for i occurrences of d (for i <= min(j, level)).
    Computed with precision p (= IndexToPrec[m]), in the same order
    as in the SageMath code.
    """
    t = mpf_mul(touslescoeffs[m][j], _R(betas[0][m], p), p, _RND)
    if j >= 1:
        t = mpf_add(t, mpf_mul(touslescoeffs[m][j-1], _R(betas[1][m], p),
                               p, _RND), p, _RND)
    if j >= 2:
        t = mpf_add(t, mpf_mul(touslescoeffs[m][j-2], _R(betas[2][m], p),
                               p, _RND), p, _RND)
    if (level > 2) and (j >= 3):
        t = mpf_add(t, mpf_mul(touslescoeffs[m][j-3], _R(betas[3][m], p),
                               p, _RND), p, _RND)
    if (level == 4) and (j >= 4):
        t = mpf_add(t, mpf_mul(touslescoeffs[m][j-4], _R(betas[4][m], p),
                               p, _RND), p, _RND)
    return t


@_fillin_irwin_docstring()
def irwin(b, d, k,
          nbdigits=34,
          level=3,
          PrecStep=500,
          all=False,
          showtimes=False,
          verbose=False,
          persistentpara=True,
          Mmax=-1
          ):
    """Somme d'Irwin pour b, d, k avec nbdigits chiffres décimaux (en tout).

    Utilise l'algorithme de Burnol, série alternée de niveau 2, 3 ou 4.

    {0}
    """

    assert 1 < level <= 4, "Le niveau (level) doit être 2 ou 3 ou 4"

    assert b > 1, "%s doit être au moins 2" % b
    bmoinsun = b - 1

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    workers = _get_workers()

    if showtimes:
        print("Préparation des précisions...", end=" ", flush=True)
        starttime = time.perf_counter()

    (nbbits,
     nbbits_final,
     Mmax,
     IndexToPrec,
     NbOfPrec) = _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax)
    P = nbbits  # precision of Rmax

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} précision(s) de précision maximale {nbbits},")
        print(f"décrémentée par multiples de {PrecStep}")

    if showtimes:
        print("Calcul des gammas...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    # lesgammas[j] is only ever needed to compute a u_{k;m} for m
    # at least equal to j. It is used only in a sum with
    # non-negative contributions, there are no subtractions. So we
    # only need it to the precision needed for u_{k;m}
    # itself.
    A1 = list(range(1, b))
    if d != 0:
        A1.remove(d)
    lesgammas = _v5_gammas(A1, bmoinsun, Mmax, IndexToPrec)

    # Those are only needed for k>O.  Same remark as for
    # lesgammas[j] relative to the precision to use.
    if k > 0:
        lespuissancesded = _v5_powers(d, Mmax, IndexToPrec)
    else:
        lespuissancesded = None

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

    workers.broadcast(("setup", IndexToPrec, A1, d, k))

    if showtimes:
        if k == 0:
            print(f"Calcul des u_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des u_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.perf_counter()

    # Recursive computation of the u_{k;m}'s.
    # touslescoeffs = [[u_{0,0}, u_{1,0}, ..., u_{k,0}],
    #                  [u_{0,1}, u_{1,1}, ..., u_{k,1}],
    #                  ...
    #                  ]
    Rb = from_int(b, P, _RND)
    Rden = from_int(b * b - bmoinsun, P, _RND)
    touslescoeffs = [ [Rb] * (k+1) ]
    c1 = [ mpf_div(mpf_mul(lesgammas[1], Rb, P, _RND), Rden, P, _RND) ]
    for j in range(1, k+1):
        c1.append(mpf_div(mpf_add(mpf_mul(_plus(lesgammas[1], d, P), Rb,
                                          P, _RND),
                                  c1[-1], P, _RND),
                          Rden, P, _RND))
    touslescoeffs.append(c1)

    PascalRows = [ [1,1] ]
    useparallel = False
    _v5_para_recurrence = _v5_setup_para_recurrence(touslescoeffs,
                                                    lesgammas,
                                                    lespuissancesded,
                                                    PascalRows,
                                                    IndexToPrec,
                                                    b, bmoinsun,
                                                    k,
                                                    showtimes,
                                                    persistentpara,
                                                    False,
                                                    workers)

    m = 1
    # We have initialized touslescoeffs[0] and touslescoeffs[1]
    # We now need for m from 2 to Mmax inclusive.
    Q, R = divmod(Mmax - 1, workers.n)
    for _P in range(Q):
        useparallel = _v5_para_recurrence(m, workers.n, useparallel)
        m += workers.n
    if R > 0:
        _ = _v5_para_recurrence(m, R, useparallel)
    del PascalRows[:]

    # the workers do not need their copies of the coefficients anymore
    workers.broadcast(("setup", IndexToPrec, A1, d, k))

    if showtimes:
        stoptime = time.perf_counter()
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} "
              + f"Fini! En tout : {stoptime-starttime:.3f}s")

    if showtimes:
        print(f"Calcul des blocs d'entiers...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level > 2:
        block3 = blocks[2]
    if level > 3:
        block4 = blocks[3]

    # The integers with level digits, according to their d-counts.
    maxblock = blocks[-1]

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={workers.n} ...")
        _lesbetas_par_nb_occurrences = _v5_map_beta_withtimes(Mmax,
                                                              workers,
                                                              maxblock)
    else:
        _lesbetas_par_nb_occurrences = _v5_map_beta_notimes(Mmax,
                                                            workers,
                                                            maxblock)

    # According to Theorem 1, formula (1) of arXiv:2402.09083, to
    # compute the m th term of the Burnol series for the Irwin sum
    # associated to exactly j occurrences we need to combine u_{j;m},
    # u_{j-1;m}, u_{j-2;m}, ... with weights which are the sum of the
    # inverse (m+1)-powers of the integers with level digits having
    # respectively 0, 1, 2, ... occurrences of digit d.
    # lesbetas[i][m] is the sum of the 1/n**(m+1) where n has level
    # digits and exactly i of them are equal to d.
    lesbetas = [ _lesbetas_par_nb_occurrences(0) ]
    if k >= 1:
        lesbetas.append(_lesbetas_par_nb_occurrences(1))
    if k >= 2:
        lesbetas.append(_lesbetas_par_nb_occurrences(2))
    if (k >= 3) and (level > 2):
        lesbetas.append(_lesbetas_par_nb_occurrences(3))
    if (k >= 4) and (level > 3):
        lesbetas.append(_lesbetas_par_nb_occurrences(4))

    workers.broadcast(("clear",))

    # Boucle qui évalue également, si all=True la série pour les j<k.
    Sk = []

    for j in range(0 if all else k, k+1):
        if showtimes:
            print(f"Calcul de l'approximation principale avec k={j}...",
                  end = ' ', flush = True)
            starttime = time.perf_counter()

        # Calcul de la série alternée de Burnol.

        S = 0

        # The Burnol formula starts with the harmonic sum of all
        # positive integers having strictly less than level digits
        # AND exactly j occurrences of the digit d.  We start with
        # the length-1 integers.  So they contribute only if j is 0
        # or 1.
        if j == 0:
            S = _suminv(block1[0], P)
        elif j == 1:
            if d != 0:
                S = mpf_div(fone, from_int(d, P, _RND), P, _RND)

        if verbose:
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(_fmt(S, P))

        # We add contribution of the length-2 integers.  This regards only
        # j=0, 1, or 2 and level must be >2.
        if 2 < level:
            if j <= 2:
                S = _plus(S, _suminv(block2[j], P), P)
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(_fmt(S, P))

        # We add contribution of the length-3 integers.  This regards only
        # j=0, 1, 2 or 3 and level must be >3.
        if 3 < level:
            if j <= 3:
                S = _plus(S, _suminv(block3[j], P), P)
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(_fmt(S, P))

        # The next contribution in the Burnol series is b times the
        # harmonic sum of length=level integers having *at most* j
        # occurrences of digit d.
        H = 0
        for i in range(1 + min(j,level)):
            H = _plus(H, _suminv(maxblock[i], P), P)
        S = _plus(S, _times(b, H, P), P)

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(_fmt(S, P))

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série alternée")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.perf_counter()

        # We now compute the Burnol series which is the alternating series
        # as in equation (1) (Theorem 1) of arXiv:2402.09083.
        # Each integer n with level digits and i occurrences of d
        # contributes u_{j-i;m}/n**(m+1).

        # We start with the smallest term contributing to the series
        # Its sign will be set later.
        bubu = _v5_series_term(touslescoeffs, lesbetas, Mmax, j, level,
                               IndexToPrec[Mmax])

        if verbose:
            lastterm = mpf_neg(bubu) if Mmax&1 else bubu
            if to_float(lastterm) == 0.:
                u, E = _v5_shorten_small_real(lastterm)
                print("The %sth term is about %f times 10^%s and it represents" % (Mmax, u, E))
            else:
                print("The %sth term is about %.3e and it represents" % (Mmax, to_float(lastterm)))

        # COMPUTATION OF THE MAIN SERIES BUILDING UP FROM SMALLEST TERMS
        # p will be the precision. When m decreases p changes from time to
        # time regularly and automatically to use more bits.
        for m in range(Mmax-1, 0, -1):  # last one is m=1
            p = IndexToPrec[m]
            # Extend the partial sum obtained to a higher precision
            # to prepare addition of a new term, which is computed
            # with higher precision.
            bubu = mpf_add(mpf_neg(bubu),
                           _v5_series_term(touslescoeffs, lesbetas, m, j,
                                           level, p),
                           p, _RND)

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = mpf_add(_R(S, P), mpf_neg(bubu), P, _RND)

        if verbose:
            pr = min(P, IndexToPrec[Mmax])
            ratio = mpf_div(lastterm, S, pr, _RND)
            if to_float(ratio) == 0.:
                u, E = _v5_shorten_small_real(ratio)
                print("%.3f 10^%s of the total." % (u, E))
            else:
                print("%.3e of the total." % to_float(ratio))

            print("La somme de m=1 à %s vaut" % Mmax)
            print(_fmt(mpf_neg(bubu), IndexToPrec[1]))

        if all:
            Sk.append(S)

    if all:
        for j in range(k+1):
            print(f"(k={j}) {_make_result(Sk[j], nbbits_final)}")

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    return _make_result(S, nbbits_final)


@_fillin_irwin_docstring()
def irwinpos(b, d, k,
             nbdigits=34,
             level=3,
             PrecStep=500,
             all=False,
             showtimes=False,
             verbose=False,
             persistentpara=True,
             Mmax=-1
             ):
    """Somme d'Irwin pour b, d, k avec nbdigits chiffres décimaux (en tout).

    Utilise algorithme de Burnol, série positive de niveau 2, 3 ou 4.

    {0}
    """

    assert 1 < level <= 4, "Le niveau (level) doit être 2 ou 3 ou 4"

    assert b > 1, "%s doit être au moins 2" % b
    bmoinsun = b - 1

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    workers = _get_workers()

    if showtimes:
        print("Préparation des précisions...", end=" ", flush=True)
        starttime = time.perf_counter()

    (nbbits,
     nbbits_final,
     Mmax,
     IndexToPrec,
     NbOfPrec) = _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax)
    P = nbbits  # precision of Rmax

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} précision(s) de précision maximale {nbbits},")
        print(f"décrémentée par multiples de {PrecStep}")

    if showtimes:
        print("Calcul des gammas ...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    # ATTENTION que la série positive a des récurrences avec b-1-d
    # à la place de d
    dprime = b - 1 - d
    A1prime = list(range(1, b))
    if dprime != 0:
        A1prime.remove(dprime)
    lesgammasprime = _v5_gammas(A1prime, bmoinsun, Mmax, IndexToPrec)

    if k > 0:
        lespuissancesdedprime = _v5_powers(dprime, Mmax, IndexToPrec)
    else:
        lespuissancesdedprime = None

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

    workers.broadcast(("setup", IndexToPrec, A1prime, dprime, k))

    if showtimes:
        if k == 0:
            print(f"Calcul des v_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des v_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.perf_counter()

    # Recursive computation of the v_{k;m}'s.  See comments in irwin().
    Rb = from_int(b, P, _RND)
    Rden = from_int(b * b - bmoinsun, P, _RND)
    touslescoeffs = [ [Rb] * (k+1) ]
    # ATTENTION: this b * b  extra is needed for the  v_{0;1}.
    c1 = [ mpf_div(_plus(b * b, mpf_mul(lesgammasprime[1], Rb, P, _RND), P),
                   Rden, P, _RND) ]
    for j in range(1, k+1):
        c1.append(mpf_div(mpf_add(mpf_mul(_plus(lesgammasprime[1], dprime, P),
                                          Rb, P, _RND),
                                  c1[-1], P, _RND),
                          Rden, P, _RND))
    touslescoeffs.append(c1)

    PascalRows = [ [1,1] ]
    useparallel = False
    _v5_para_recurrence = _v5_setup_para_recurrence(touslescoeffs,
                                                    lesgammasprime,
                                                    lespuissancesdedprime,
                                                    PascalRows,
                                                    IndexToPrec,
                                                    b, bmoinsun,
                                                    k,
                                                    showtimes,
                                                    persistentpara,
                                                    True,
                                                    workers)
    m = 1
    # We have initialized touslescoeffs[0] and touslescoeffs[1]
    # We now need for m from 2 to Mmax inclusive.
    Q, R = divmod(Mmax - 1, workers.n)
    for _P in range(Q):
        useparallel = _v5_para_recurrence(m, workers.n, useparallel)
        m += workers.n
    if R > 0:
        _= _v5_para_recurrence(m, R, useparallel)
    del PascalRows[:]

    # the workers do not need their copies of the coefficients anymore
    workers.broadcast(("setup", IndexToPrec, A1prime, dprime, k))

    if showtimes:
        stoptime = time.perf_counter()
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} (fait) "
              + f"{stoptime-starttime:.3f}s")

    if showtimes:
        print(f"Calcul des blocs d'entiers...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level > 2:
        block3 = blocks[2]
    if level > 3:
        block4 = blocks[3]
    maxblock = blocks[-1]
    # ATTENTION!
    # We shift by +1 all integers in sublists of maxblock. This is
    # to avoid having to use n+1 afterwards for inverse power sums.
    maxblockshifted = [[ n + 1  for n in L] for L in maxblock]

    # Pay attention that _lesbetas_par_nb_occurrences produces here
    # beta's which are sums of 1/(n+1)**(m+1)'s for certain n's
    # whereas in irwin() it was sums of 1/n**(m+1).
    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={workers.n} ...")
        _lesbetas_par_nb_occurrences = _v5_map_beta_withtimes(Mmax,
                                                              workers,
                                                              maxblockshifted)
    else:
        _lesbetas_par_nb_occurrences = _v5_map_beta_notimes(Mmax,
                                                            workers,
                                                            maxblockshifted)

    # lesbetas[i][m] is the sum of the 1/(n+1)**(m+1) where n has level
    # digits and exactly i of them are equal to d.
    lesbetas = [ _lesbetas_par_nb_occurrences(0) ]
    if k >= 1:
        lesbetas.append(_lesbetas_par_nb_occurrences(1))
    if k >= 2:
        lesbetas.append(_lesbetas_par_nb_occurrences(2))
    if (k >= 3) and (level > 2):
        lesbetas.append(_lesbetas_par_nb_occurrences(3))
    if (k >= 4) and (level > 3):
        lesbetas.append(_lesbetas_par_nb_occurrences(4))

    workers.broadcast(("clear",))

    # Boucle qui évalue également, si all=True la série pour les j<k.
    Sk = []

    for j in range(0 if all else k, k+1):
        if showtimes:
            print(f"Calcul de l'approximation principale avec k={j}...",
                  end = ' ', flush = True)
            starttime = time.perf_counter()

        # Calcul de la série positive de Burnol.
        S = 0

        if j == 0:
            S = _suminv(block1[0], P)
        elif j == 1:
            if d != 0:
                S = mpf_div(fone, from_int(d, P, _RND), P, _RND)

        if verbose:
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(_fmt(S, P))

        if 2 < level:
            if j <= 2:
                S = _plus(S, _suminv(block2[j], P), P)
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(_fmt(S, P))

        if 3 < level:
            if j <= 3:
                S = _plus(S, _suminv(block3[j], P), P)
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(_fmt(S, P))

        # ATTENTION
        #
        # So far the S value is the same as for irwin(). But for the
        # positive series, the next contribution in the Burnol
        # series is b times the sum of the 1/(n+1) (not 1/n) where
        # the n's have exactly "level" digits and *at most* j
        # occurrences of digit d.  This is why we have the
        # maxblockshifted[i] here.
        H = 0
        for i in range(1 + min(j,level)):
            H = _plus(H, _suminv(maxblockshifted[i], P), P)
        S = _plus(S, _times(b, H, P), P)

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(_fmt(S, P))

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série positive")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.perf_counter()

        # We now compute the Burnol series which is the positive series
        # as in equation (4) (Theorem 4) of arXiv:2402.09083.
        # Each integer n with level digits and i occurrences of d
        # contributes v_{j-i;m}/(n+1)**(m+1).
        # We start with the smallest term contributing to the series
        bubu = _v5_series_term(touslescoeffs, lesbetas, Mmax, j, level,
                               IndexToPrec[Mmax])

        if verbose:
            lastterm = bubu
            if to_float(lastterm) == 0.:
                u, E = _v5_shorten_small_real(lastterm)
                print(f"The {Mmax}th term is about {u:f} 10^{E} i.e. ",
                      end = "", flush = True)
            else:
                print(f"The {Mmax}th term is about {to_float(lastterm):.3e} i.e. ",
                      end = "", flush = True)

        # COMPUTATION OF THE MAIN SERIES BUILDING UP FROM SMALLEST TERMS
        for m in range(Mmax-1, 0, -1):  # last one is m=1
            p = IndexToPrec[m]
            bubu = mpf_add(bubu,
                           _v5_series_term(touslescoeffs, lesbetas, m, j,
                                           level, p),
                           p, _RND)

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = mpf_add(_R(S, P), bubu, P, _RND)

        if verbose:
            pr = min(P, IndexToPrec[Mmax])
            ratio = mpf_div(lastterm, S, pr, _RND)
            if to_float(ratio) == 0.:
                u, E = _v5_shorten_small_real(ratio)
                print("%.3f 10^%s of the total." % (u, E))
            else:
                print("%.3e of the total." % to_float(ratio))

            print("La somme de m=1 à %s vaut" % Mmax)
            print(_fmt(bubu, IndexToPrec[1]))

        if all:
            Sk.append(S)

    if all:
        for j in range(k+1):
            print(f"(k={j}) {_make_result(Sk[j], nbbits_final)}")

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    return _make_result(S, nbbits_final)


if __name__ == "__main__":
    print(f"""
Hello, this file {__filename__} provides two functions irwin()
and irwinpos().  Use help(irwin) or help(irwinpos) for help.

This is version {__version__} of {__date__}.

General information is also available in the irwin_v5_docstring
variable: print(irwin_v5_docstring).

The variable "maxworkers" sets the number of worker processes.
You can set it via irwin_v5_mpmath_spawn.maxworkers = N
(no need to reload).

{maxworkersinfostring}
"""
          )
