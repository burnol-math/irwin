# -*- mode: python ; coding: utf-8; -*-

# irwin_v5.sage
# Use via load("irwin_v5.sage") in sage interactive mode

__version__  = "1.5.7"
__date__     = "2025/05/17"
__filename__ = "irwin_v5.sage"

irwin_v5_docstring = """
This file is an evolution of irwin.sage as available at

    https://arxiv.org/src/2402.09083v5/anc

with three main enhancements relevant to targeting 1000 or
more decimal digits of precision:

- Use of decreasing precision for higher terms in the series.
  This is done with a granularity of PrecStep bits.  It is an
  optional parameter for the procedures irwin() and irwinpos().
  It defaults to 500.

- There is no pre-computation of the first 1000 rows of the Pascal
  triangle anymore.  Only maxworkers (see item on parallelization
  next) rows of the Pascal triangle are kept at any given time in
  memory.  See file taille_pascal.pdf for details on the storage
  size needed for rows of the Pascal triangle.

- Parallelization, using @parallel(ncpus=maxworkers), where
  maxworkers defaults to 8 and can be defined prior to loading
  this file (no attempt is made to adjust dynamically to the
  number of cores available).  Users of macOS are advised to
  look at issue #1 at the repository and check if it applies
  to their system.

Miscellaneous remarks:

* The default level is now 3 (it used to be 2). When one goes
  into hundreds of digits it is advantageous, and for less
  digits the computation is fast enough anyhow (info: for circa
  less than 100-150 decimal digits one can get faster execution
  times from using the irwin_v3.sage provided irwin(), which
  uses level=2 per default and has no parallelization).

* The original 2024 version, if used with "all = True" printed
  with 2 or 3 extra decimal digits the intermediate Irwin sums.
  Now, all printed values round away the nbguardbits guard bits
  used internally for the computations.

* The global variable nbguardbits initial value was increased
  from 8 to 12.

* The way the number of terms of the series to keep is chosen was
  slightly modified and usually one more term is kept (for
  b=10). The new way is safer for small b's or higher k's.

Thanks to Arnaud Bodin and Yusuf Emin Akpinar for interactions.

Copyright (C) 2025 Jean-François Burnol
License: CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/
"""

# https://stackoverflow.com/a/10308363
# def _docstring_parameter(*sub):
#     def dec(obj):
#         obj.__doc__ = obj.__doc__.format(*sub)
#         return obj
#     return dec

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
        Terms of the series are computed with a RealField of
        evolving precision, which differs from the maximal
        precision by a suitable multiple of PrecStep.
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

    :rtype: :class:`sage.rings.real_mpfr.RealNumber`
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b.

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
    set maxworkers appropriately prior to (re-) loading the
    module.

    Example:
    --------

    sage: {1}(10, 9, 4, 52, all=True)
    (k=0) 22.92067661926415034816365709437593191494476243699848
    (k=1) 23.04428708074784831967594930973617482538959203064774
    (k=2) 23.02604026596124378845022249787272342108112267542086
    (k=3) 23.02585299837244431714290384468012275518705238435290
    (k=4) 23.02585095265829261377053973815542996035002267989413
    23.02585095265829261377053973815542996035002267989413
"""

import time
nbguardbits = 12

try:
    _ = maxworkers
    maxworkersinfostring = ("maxworkers variable already existed "
                            f"and has value {maxworkers}.")
except NameError:
    maxworkers = 8
    maxworkersinfostring = ("maxworkers variable has been created "
                            "and assigned value 8.")


@parallel(ncpus=maxworkers)
def _v5_ukm_partial(a, m, Pm, G, D, T, Rm, k):
    """Recurrences (partial) for the u_{j;m}'s or v_{j;m}'s.

    - This handles all j's from 0 to k (because to compute
      for k we need for k-1, and to compute for k-1, we
      need for k-2 and so on until j=0).
    - Pm[i] stands for binomial coefficient "m choose i".
    - G[i] stands for gamma (or gammaprime) power sum.
    - D[i] is d**i or dprime**i (dprime = b-1-d)
    - T[n][j] holds previously known u_{j;n}'s or v_{j;n}'s.
    - Rm is a RealField using a precision suitable for the
      evaluation of the u_{j;m}'s or v_{j;m}'s.

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
    A = list(sum(Pm[i]*Rm(G[i])*Rm(T[m - i][j]) for i in range(a, m + 1))
             for j in range(k + 1))
    B = [ 0 ]
    B.extend(sum(Pm[i]*Rm(D[i])*Rm(T[m - i][j-1]) for i in range(a, m + 1))
             for j in range(1, k + 1))
    return [ A[j] + B[j] for j in range(k + 1) ]


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
                              PascalRows, IndexToR, b, bmoinsun, k,
                              showtimes, persistentpara, is_for_vm):
    """Set up procedure calling _v5_ukm_partial and completing its job.
    """
    def _v5_para_recurrence(m, step, useparallel):
        """Wrapper of parallelized calls to _v5_ukm_partial().

        First we compute "step" (which is maxworkers or less than it)
        new rows of the Pascal triangle of binomial coefficients.  We
        use some specificities of how Python handles list type to do
        that in a way persistent in memory across calls.

        Then, if useparallel is True we call the parallelized
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
        if ((M - 400) % 500 < maxworkers):
            starttime_ns = time.perf_counter_ns()
            results = _v5_ukm_partial(((a, M + a,
                                        PascalRows[a],
                                        Gammas,
                                        PuissancesDeD,
                                        touslescoeffs,
                                        IndexToR[M + a],
                                        k)
                                       for a in range(1, step + 1)))
            ukm_partial = [ None ]
            ukm_partial.extend([result[1] for result in sorted(list(results))])
            multitime_ns = time.perf_counter_ns() - starttime_ns

            if useparallel and persistentpara:
                # do not check again
                if showtimes:
                    print(f"... mode parallèle persistant ({multitime_ns*1e-9:.3f}s; "
                          f"{M}<m<={M+step})")
            else:
                starttime_ns = time.perf_counter_ns()
                m = M
                ukm_partial = [ None ]
                for j in range(1, 1 + step):
                    m += 1
                    Rm = IndexToR[m]
                    Am = [ sum(PascalRows[j][i]
                               * Rm(Gammas[i])
                               * Rm(touslescoeffs[m-i][p])
                               for i in range(j, m+1)) for p in range(k+1) ]
                    Bm = [0]
                    Bm.extend([ sum(PascalRows[j][i]
                                    * Rm(PuissancesDeD[i])
                                    * Rm(touslescoeffs[m-i][p-1])
                                    for i in range(j, m+1)) for p in range(1, k+1) ])
                    ukm_partial.append([ Am[n] + Bm[n] for n in range(k+1) ])
                singletime_ns = time.perf_counter_ns() - starttime_ns

                if showtimes:
                    _v5_umtimeinfo(singletime_ns, multitime_ns,
                                   useparallel,
                                   maxworkers, M, step)

                useparallel = multitime_ns < singletime_ns

        elif useparallel:
            results = _v5_ukm_partial(((a, M + a,
                                        PascalRows[a],
                                        Gammas,
                                        PuissancesDeD,
                                        touslescoeffs,
                                        IndexToR[M + a],
                                        k)
                                       for a in range(1, step + 1)))
            ukm_partial = [ None ]
            ukm_partial.extend([result[1] for result in sorted(list(results))])

        else:
            m = M
            ukm_partial = [ None ]
            for j in range(1, 1 + step):
                m += 1
                Rm = IndexToR[m]
                Am = [ sum(PascalRows[j][i]
                           * Rm(Gammas[i])
                           * Rm(touslescoeffs[m-i][p])
                           for i in range(j, m+1)) for p in range(k+1) ]
                Bm = [0]
                Bm.extend([ sum(PascalRows[j][i]
                                * Rm(PuissancesDeD[i])
                                * Rm(touslescoeffs[m-i][p-1])
                                for i in range(j, m+1)) for p in range(1, k+1) ])
                ukm_partial.append([ Am[n] + Bm[n] for n in range(k+1) ])

        # Now correct the um's (or vm's) (prior to dividing by b**(m+1)-b+1)
        # via the addition of finitely missing contributions in order of increasing
        # m's.
        m = M
        for j in range(1, 1 + step):
            m += 1
            Rm = IndexToR[m]
            D = Rm( b**(m+1) - bmoinsun )
            # Attention to the b**(m+1) extra term specific to v_m recurrence.
            # Attension that parentheses are needed to delimit what "else" caches.
            cm = [ ((Rm(b ** (m+1)) if is_for_vm else 0)
                    + ukm_partial[j][0]
                    + sum(PascalRows[j][i]
                          * Rm(Gammas[i])
                          * Rm(touslescoeffs[m-i][0])
                          for i in range(1, j))
                    ) / D
                  ]
            for p in range(1, k+1):
                _ = (ukm_partial[j][p]
                     + sum(PascalRows[j][i]
                           * Rm(Gammas[i])
                           * Rm(touslescoeffs[m-i][p])
                           for i in range(1, j))
                     + cm[-1]
                     + sum(PascalRows[j][i]
                           * Rm(PuissancesDeD[i])
                           * Rm(touslescoeffs[m-i][p-1])
                           for i in range(1, j))
                     ) / D
                cm.append(_)
            touslescoeffs.append(cm)
        # Update status.
        return useparallel
    return _v5_para_recurrence


def _v5_beta_aux(m, R, nblock):
    return sum(1/R(n ** (m+1)) for n in nblock)


@parallel(ncpus=maxworkers)
def _v5_beta(start, end, IR, nblock):
    """Parallelized caller to computation of beta coefficients.

    For m varying in a given range via steps of value maxworkers, we
    compute the sum of 1/n**(m+1) for n varying in given "nblock".
    IR stands for IndexToR which maps indices m to suitable
    RealField specifying the used precision.  Higher indices use
    lower precision, this is why the range of m's is split according
    to value modulo maxworkers, so that the computation costs are
    about equal across workers.

    The start will be an integer from 1 (not zero) to maxworkers.
    The end is simply Mmax+1, so the last index m used is Mmax.
    """
    return list(_v5_beta_aux(m, IR[m], nblock)
                for m in range(start, end, maxworkers))


def _v5_map_beta_notimes(Mmax, IndexToR, maxblock):
    """Sets up a procedure to call _v5_beta() and assembles its results.

    The defined procedure will receive an argument j which is in the
    range from 0 to k inclusive.  It will then use the integers in
    the "block" maxblock[j] as the ones for which the sum of inverse
    powers needs to be computed.

    The procedure defined by this does not display intermediate
    computing times.  Depending on whether Mmax is a multiple of
    maxworkers or not two procedures are defined, but this is a bit
    silly because the gain is minuscule as the defined procedures
    will be called only k+1 times.
    """
    extra = maxworkers - ( Mmax % maxworkers )
    if extra < maxworkers:
        def map__v5_beta(j):
            """Calls parallelized _v5_beta() and assembles its results.

            After having computed beta_{m+1}'s for m's split by
            their modulo maxworkers value (in (1,..., maxworkers))
            we reorganize the maxworkers lists of values into a
            single list in order of increasing m's.

            When Mmax is not a multiple of maxworkers, the returned
            lists have two distinct lengths, and before zipping we
            extend the shorter ones by None.  Zipping will then have
            a number of extra None's at the end which we then
            remove.  The syntax used appears to be Pythonic, I don't
            know it is the most efficient and if going via zip() is
            good idea.  As said, it is Pythonic at least.
            """
            L = [0]
            inputdata = [(i,
                          Mmax + 1,
                          IndexToR,
                          maxblock[j])
                         for i in range(1, 1 + maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v5_beta(inputdata)))]
            for j in range(1, extra + 1):
                results_1[-j].append(None)
            L.extend(x for xs in zip(*results_1) for x in xs)
            return L[:-extra]
    else:
        def map__v5_beta(j):
            """Calls parallelized _v5_beta() and assembles its results.

            After having computed beta_{m+1}'s for m's split by
            their modulo maxworkers value (in (1,..., maxworkers))
            we reorganize the maxworkers lists of values into a
            single list in order of increasing m's.
            """
            L = [0]
            inputdata = [(i,
                          Mmax + 1,
                          IndexToR,
                          maxblock[j])
                         for i in range(1, 1 + maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v5_beta(inputdata)))]
            L.extend(x for xs in zip(*results_1) for x in xs)
            return L
    return map__v5_beta


def _v5_map_beta_withtimes(Mmax, IndexToR, maxblock):
    """Sets up a procedure to call _v5_beta() and assembles its results.

    The defined procedure will receive an argument j which is in the
    range from 0 to k inclusive.  It will then use the integers in
    the "block" maxblock[j] as the ones for which the sum of inverse
    powers needs to be computed.

    The procedure defined by this does displays intermediate
    computing times.  It will divide fomr this the range from 1 to
    Mmax in chunks of size a multiple of maxworkers near to 1000.
    If maxworkers if 32 or more, chunks of size 32*maxworkers are
    used for displaying their timings.
    """
    # We want to display some visual sign of progress.
    # Find the largest multiple of maxworkers at most 1000,
    # do something reasonable if maxworkers is big
    q = max(1000 // maxworkers, 32)
    mSize = q * maxworkers
    def map__v5_beta(j):
        """Calls parallelized _v5_beta() and assembles its results.

        And compute intermediate timings while doing it.

        After having computed beta_{m+1}'s for m's split by their
        modulo maxworkers value (in (1,..., maxworkers)) in various
        ranges we need to reorganize the maxworkers lists of values
        in order of increasing m's and extend the list which will
        hold all the values.
        """
        print(f"... ({j} occ.) ", end = "", flush = True)
        starttime = time.perf_counter()
        lasttime = starttime
        mbegin = 1  # will remain congruent to 1 modulo maxworkers
        mend = 1    # this one also
        L = [0]
        for rep in range(Mmax // mSize):
            mend   = mbegin + mSize
            # In this loop, mSize is a multiple of q.
            # We call the parallelized _v5_beta with exactly maxworkers
            # arguments.
            # Memo: le premier argument décidera du sorted.
            inputdata = [(mbegin + i,
                          mend,
                          IndexToR,
                          maxblock[j])
                         for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v5_beta(inputdata)))]
            L.extend(x for xs in zip(*results_1) for x in xs)
            stoptime = time.perf_counter()
            print(f"m<{mend} ({stoptime-lasttime:.3f}s)",
                  end = "\n             ", flush= True)
            lasttime = stoptime
            mbegin = mend

        if mend < Mmax+1:
            extra = maxworkers - ( (Mmax + 1 - mend) % maxworkers )
            inputdata = [(mend + i,
                          Mmax + 1,
                          IndexToR,
                          maxblock[j])
                         for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v5_beta(inputdata)))]
            if extra > 0:
                for j in range(1, extra + 1):
                    results_1[-j].append(None)
                L.extend(x for xs in zip(*results_1) for x in xs)
                del L[-extra:]
            else:
                L.extend(x for xs in zip(*results_1) for x in xs)
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
    # I have to work around better this problem with float conversion
    # possibly giving zero due to low exponent
    x = RealField(53)(rr)
    s, _, _ = x.sign_mantissa_exponent()
    y = abs(x).log10()
    N = y.trunc()
    return s*10**float(y-N+1), N-1


def _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax=-1):
    """Preparation of an array mapping each m to a RealField.

    See irwin_v5_doc.pdf for mathematical details.
    """

    # Chose number of bits to (try to) guarantee we will have nbdigits
    # decimal digits in output.
    nbbits_final = ceil((nbdigits+1)*log(10,2))
    Rfinal = RealField(nbbits_final)

    # Computations are done (for the main terms) with elevated precision.
    nbbits = nbbits_final + nbguardbits
    R = RealField(nbbits)

    # See irwin_v5_doc.pdf for the mathematical justification for this
    # choice of Mmax, which is the number of terms used from the
    # series given in Burnol papers.
    _Mmax = floor((nbbits - nbguardbits/2)/(level-1)/log(b,2))

    # The number of distinct precisions we need.
    NbOfPrec = 1 + floor((nbbits - nbguardbits/2)/PrecStep)
    # We pre-create all needed RealField's and store them in a
    # list.
    LesReels = [R]
    for j in range(1, NbOfPrec):
        LesReels.append(RealField(nbbits - j * PrecStep))

    # See irwin_v5_doc.pdf for the justification that we only need
    #
    #    nbbits - (l-1) * m * log(b,2)
    #
    # precision for the computation of the mth term.  So j is
    # chosen to be the largest such that nbbits - jT is at least
    # that value.  Hence j is floor((l-1)*log(b,2)*m/T) (with T =
    # PrecStep).

    # To avoid having to compute this j for each m, we prepare a
    # list IndexToR which will hav pre-computed this for all used
    # m's.  This all is probably very silly.
    IndexToR = []
    #
    # We need to know which m will use a given j.
    # They verify j<= (l-1)*log(b,2)*m/T < j+1, i.e.
    # j*T/(l-1)/log(b,2)<= m < (j+1)*T/(l-1)/log(b,2), i.e.:
    #
    # ceil(j*T/(l-1)/log(b,2))<= m < ceil((j+1)*T/(l-1)/log(b,2)).
    #
    # If T/(l-1)/log(b,2)<1 it is not guaranteed that a given j will
    # have associated m's. But this is not a problem.
    IndexToPrecRatio = PrecStep / (level - 1) / log(b,2)  # not an integer !
    oldindexbound = 0
    for j in range(NbOfPrec):
        newindexbound = ceil((j+1)*IndexToPrecRatio)
        IndexToR.extend([LesReels[j]] * (newindexbound - oldindexbound))
        oldindexbound = newindexbound

    # For j=0, newindexbound is at least 1 but may be not at least 2,
    # But we want m=1 to use always maximal precision.
    IndexToR[1] = R  # LesReels[0]
    # In total the length of this array is ceil(NbOfPrec*IndexToPrecRatio).
    # It is at least
    #         ceil((nbbits - nbguardbits/2)/(level-1)/log(b,2)).
    # On the other hand
    #    Mmax=floor((nbbits - nbguardbits/2)/(level-1)/log(b,2)).
    # For the two to be the same we would need
    # (nbbits-nbguardbits/2)/(level-1)/log(b,2) to be an
    # integer. But this integer would be strictly less than
    # NbOfPrec*IndexToPrecRatio (hence less than the ceil())
    # because NbOfPrec is strictly more
    # than (nbbits - nbguardbits/2)/PrecStep, so in fact the
    # length of the array is always at least Mmax + 1 which is
    # what we want.  But as Python float computations with reals
    # are not exact, we do this additional precautionary step.
    IndexToR.extend([LesReels[-1]] * (_Mmax + 1 - len(IndexToR)))

    if Mmax == -1:
        Mmax = _Mmax
    else:
        # If user specified a custom Mmax which is beyond our
        # estimate we need to extend.
        if Mmax > _Mmax:
            IndexToR.extend([LesReels[-1]] * (Mmax - _Mmax))
            print(f"!!!! Warning Mmax={Mmax} is probably needlessly big.")
            print(f"!!!! {_Mmax} should be enough but we will use your value.")
            print(f"!!!! Use verbose=True to see the size of the smallest term.")

    return nbbits, nbbits_final, R, Rfinal, Mmax, IndexToR, NbOfPrec


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

    # SURTOUT NE PAS FAIRE if type(b) == type(1) !!!!
    # Ça marche pour des inputs directs mais pas pour "for b in range(B)".
    assert b > 1, "%s doit être au moins 2" % b
    bmoinsun = b - 1

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    if showtimes:
        print("Préparation des RealField...", end=" ", flush=True)
        starttime = time.perf_counter()

    (nbbits,
     nbbits_final,
     Rmax,
     Rfinal,
     Mmax,
     IndexToR,
     NbOfPrec) = _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax)

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} RealField(s) de précision maximale {nbbits},")
        print(f"décrémentée par multiples de {PrecStep}")

    if showtimes:
        print("Calcul des gammas...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    # lesgammas[j] is only ever needed to compute a u_{k;m} for m
    # at least equal to j. It is used only in a sum with
    # non-negative contributions, there are no subtractions. So we
    # only need it to the precision needed for u_{k;m}
    # itself. There is however a cumulative error effect in a sum
    # with m contributions, but it should be absorbed in the
    # safety cushion we have in place in evaluating the needed
    # precision.  Perhaps if k is very large, this safety cushion
    # could prove defective.
    A1 = list(range(1, b))
    if d != 0:
        A1.remove(d)
    lesgammas = [ bmoinsun ]
    for m in range(1, Mmax+1):
        Rm = IndexToR[m]
        lesgammas.append(Rm(sum(a**m for a in A1)))

    # Those are only needed for k>O.  Same remark as for
    # lesgammas[j] relative to the precision to use.
    if k > 0:
        lespuissancesded = [ 1 ]
        for m in range(1, Mmax+1):
            Rm = IndexToR[m]
            lespuissancesded.append(Rm(d**m))
    else:
        # we need the name to be defined when calling _v5_ukm_partial
        lespuissancesded = None

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

    if showtimes:
        if k == 0:
            print(f"Calcul des u_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des u_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.perf_counter()

    # Recursive computation of the u_{k;m}'s.
    # In order to have to evaluate each Pascal triangle row only once, we
    # use a slightly modified syntax compared to the 2024 version, instead
    # of having touslescoeffs = [ [ the u_{0,m}'s ], [ the u_{1,m}'s ], ... ]
    # it is now touslescoeffs = [[u_{0,0}, u_{1,0}, ..., u_{k,0}],
    #                            [u_{0,1}, u_{1,1}, ..., u_{k,1}],
    #                            ...
    #                            ]
    touslescoeffs = [ [Rmax(b)] * (k+1) ]
    c1 = [ lesgammas[1] * Rmax(b) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(((lesgammas[1] + d) * Rmax(b) + c1[-1])/Rmax(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRows = [ [1,1] ]
    useparallel = False
    _v5_para_recurrence = _v5_setup_para_recurrence(touslescoeffs,
                                                    lesgammas,
                                                    lespuissancesded,
                                                    PascalRows,
                                                    IndexToR,
                                                    b, bmoinsun,
                                                    k,
                                                    showtimes,
                                                    persistentpara,
                                                    False)

    m = 1
    # We have initialized touslescoeffs[0] and touslescoeffs[1]
    # We now need for m from 2 to Mmax inclusive.
    Q, R = divmod(Mmax - 1, maxworkers)
    for P in range(Q):
        useparallel = _v5_para_recurrence(m, maxworkers, useparallel)
        m += maxworkers
    # Ici on va invoquer une procédure parallélisée avec < maxworkers.
    if R > 0:
        _ = _v5_para_recurrence(m, R, useparallel)

    if showtimes:
        stoptime = time.perf_counter()
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} "
              + f"Fini! En tout : {stoptime-starttime:.3f}s")

    if showtimes:
        print(f"Calcul des blocs d'entiers...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    # Calcul des blocs d'entiers suivant longueur et nombre d'occurrences.
    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]  # list [ [non-zero digits not d], [d] or []]
    block2 = blocks[1]  # blocks2[j] = integers with 2 digits and j among
                        # them are equal to d.
    if level > 2:
        block3 = blocks[2]  # integers with 3 digits, assembled according
                            # to d-count.
    if level > 3:
        block4 = blocks[3]  # integers with 4 digits, according to d-count.

    # The integers with level digits, according to their d-counts.
    maxblock = blocks[-1]
    # NOTA BENE: maxblock will have as last element an empty [] if d=0
    #            This empty [] will not cause problems for the sum()'s
    #            such as sum(1/Rmax(x) for x in maxblock[i])

    # Calcul parallèle des beta (sommes d'inverses de puissances).
    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
        _lesbetas_par_nb_occurrences = _v5_map_beta_withtimes(Mmax,
                                                              IndexToR,
                                                              maxblock)
    else:
        _lesbetas_par_nb_occurrences = _v5_map_beta_notimes(Mmax,
                                                            IndexToR,
                                                            maxblock)

    # According to Theorem 1, formula (1) of arXiv:2402.09083, to
    # compute the m th term of the Burnol series for the Irwin sum
    # associated to exactly j occurrences we need to combine u_{j;m},
    # u_{j-1;m}, u_{j-2;m}, ... with weights which are the sum of the
    # inverse (m+1)-powers of the integers with level digits having
    # respectively 0, 1, 2, ... occurrences of digit d.

    # This returns the list L0 such that L0[m] is the sum of the 1/n**(m+1)
    # where n has level digits and none of them is d.
    # If showtimes is True it prints timings.
    lesbetas_maxblock0 = _lesbetas_par_nb_occurrences(0)

    # This returns the list L1 such that L1[m] is the sum of the 1/n**(m+1)
    # where n has level digits and exactly one of them is d.
    if k >= 1:
        lesbetas_maxblock1 = _lesbetas_par_nb_occurrences(1)

    # This returns the list L2 such that L2[m] is the sum of the 1/n**(m+1)
    # where n has level digits and exactly two of them are d.
    if k >= 2:
        lesbetas_maxblock2 = _lesbetas_par_nb_occurrences(2)

    # idem for 3 occurrences
    if (k >= 3) and (level > 2):
        lesbetas_maxblock3 = _lesbetas_par_nb_occurrences(3)
    # idem for 4 occurrences
    if (k >= 4) and (level > 3):
        lesbetas_maxblock4 = _lesbetas_par_nb_occurrences(4)

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
            S = sum(1/Rmax(x) for x in block1[0])
        elif j == 1:
            if d != 0:
                S = 1/Rmax(d)

        if verbose:
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(S)

        # We add contribution of the length-2 integers.  This regards only
        # j=0, 1, or 2 and level must be >2.
        if 2 < level:
            if j <= 2:
                S += sum(1/Rmax(x) for x in block2[j])
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(S)

        # We add contribution of the length-3 integers.  This regards only
        # j=0, 1, 2 or 3 and level must be >3.
        if 3 < level:
            if j <= 3:
                S += sum(1/Rmax(x) for x in block3[j])
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(S)

        # The next contribution in the Burnol series is b times the
        # harmonic sum of length=level integers having *at most* j
        # occurrences of digit d.
        # Jusqu'à j répétitions ; si j >= level, on s'arrête à
        # level répétitions max
        S += b * (sum(sum(1/Rmax(x) for x in maxblock[i])
                      for i in range(1 + min(j,level))))

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(S)

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série alternée")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.perf_counter()

        # We now compute the Burnol series which is the alternating series
        # as in equation (1) (Theorem 1) of arXiv:2402.09083.

        # We start with the smallest term contributing to the series
        # Its sign will be set later.
        Rm = IndexToR[-1]  # RealField at lowest used precision.

        # Compared to 2024 version touslescoeffs has its two indices permuted
        # Each integer n with level digits and no occurrence of d contributes
        # u_{j;m}/n**(m+1).
        bubu = touslescoeffs[Mmax][j] * lesbetas_maxblock0[Mmax]
        # Each integer n with level digits and 1 occurrence of d contributes
        # u_{j-1;m}/n**(m+1).
        if j >= 1:
            bubu += touslescoeffs[Mmax][j-1] * lesbetas_maxblock1[Mmax]
        # Each integer n with level digits and 2 occurrences of d contributes
        # u_{j-2;m}/n**(m+1).
        if j >= 2:
            bubu += touslescoeffs[Mmax][j-2] * lesbetas_maxblock2[Mmax]
        # Each integer n with level digits and 3 occurrences of d contributes
        # u_{j-3;m}/n**(m+1).
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[Mmax][j-3] * lesbetas_maxblock3[Mmax]
        # Each integer n with level digits and 4 occurrences of d contributes
        # u_{j-4;m}/n**(m+1).  This can happen here only if level is 4,
        # as we have not implemented higher levels.  And if d is non zero,
        # there is one contribution, if d is zero there are none.
        if (level == 4) and (j >= 4):
            bubu += touslescoeffs[Mmax][j-4] * lesbetas_maxblock4[Mmax]

        if verbose:
            lastterm = -bubu if Mmax&1 else bubu
            if float(lastterm) == 0.:
                u, E = _v5_shorten_small_real(lastterm)
                print("The %sth term is about %f times 10^%s and it represents" % (Mmax, u, E))
            else:
                print("The %sth term is about %.3e and it represents" % (Mmax, lastterm))

        # COMPUTATION OF THE MAIN SERIES BUILDING UP FROM SMALLEST TERMS
        # Rm will be the RealField. When m decreases Rm changes from time to
        # time regularly and automatically to use more bits.
        for m in range(Mmax-1, 0, -1):  # last one is m=1
            Rm = IndexToR[m]
            # Extend the partial sum obtained to a RealField using
            # higher precision to prepare addition of a new term,
            # which is computed with higher precision.
            bubu = Rm(-bubu)
            #
            # See comments above about the contributions 1/n**(m+1)
            # for integers n having level digits, depending on the count
            # of d's.
            bubu += touslescoeffs[m][j] * lesbetas_maxblock0[m]

            if j >= 1:
                bubu += touslescoeffs[m][j-1] * lesbetas_maxblock1[m]
            if j >= 2:
                bubu += touslescoeffs[m][j-2] * lesbetas_maxblock2[m]
            if (level > 2) and (j >= 3):
                bubu += touslescoeffs[m][j-3] * lesbetas_maxblock3[m]
            if (level == 4) and (j >= 4):
                bubu += touslescoeffs[m][j-4] * lesbetas_maxblock4[m]

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = S - bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = _v5_shorten_small_real(ratio)
                print("%.3f 10^%s of the total." % (u, E))
            else:
                print("%.3e of the total." % ratio)

            print("La somme de m=1 à %s vaut" % Mmax)
            print(-bubu)

        if all:
            Sk.append(S)

    if all:
        for j in range(k+1):
            print(f"(k={j}) {Rfinal(Sk[j])}")

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    return Rfinal(S)


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

    if showtimes:
        print("Préparation des RealField...", end=" ", flush=True)
        starttime = time.perf_counter()

    (nbbits,
     nbbits_final,
     Rmax,
     Rfinal,
     Mmax,
     IndexToR,
     NbOfPrec) = _v5_setup_realfields(nbdigits, PrecStep, b, level, Mmax)

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} RealField(s) de précision maximale {nbbits},")
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
    lesgammasprime = [ bmoinsun ]
    for m in range(1, Mmax+1):
        Rm = IndexToR[m]
        lesgammasprime.append(Rm(sum(a**m for a in A1prime)))

    if k > 0:
        lespuissancesdedprime = [ 1 ]
        for j in range(1, Mmax+1):
            Rj = IndexToR[j]
            lespuissancesdedprime.append(Rj(dprime**j))
    else:
        lespuissancesdedprime = None

    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

    if showtimes:
        if k == 0:
            print(f"Calcul des v_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des v_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.perf_counter()

    # Recursive computation of the v_{k;m}'s.  See comments in irwin().
    touslescoeffs = [ [Rmax(b)] * (k+1) ]
    # ATTENTION: this b * b  extra is needed for the  v_{0;1}.
    c1 = [ (b * b + lesgammasprime[1] * Rmax(b)) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(( (lesgammasprime[1] + dprime) * Rmax(b)
                    + c1[-1])/Rmax(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRows = [ [1,1] ]
    useparallel = False
    _v5_para_recurrence = _v5_setup_para_recurrence(touslescoeffs,
                                                    lesgammasprime,
                                                    lespuissancesdedprime,
                                                    PascalRows,
                                                    IndexToR,
                                                    b, bmoinsun,
                                                    k,
                                                    showtimes,
                                                    persistentpara,
                                                    True)
    m = 1
    # We have initialized touslescoeffs[0] and touslescoeffs[1]
    # We now need for m from 2 to Mmax inclusive.
    Q, R = divmod(Mmax - 1, maxworkers)
    for P in range(Q):
        useparallel = _v5_para_recurrence(m, maxworkers, useparallel)
        m += maxworkers
    if R > 0:
        _= _v5_para_recurrence(m, R, useparallel)

    if showtimes:
        stoptime = time.perf_counter()
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} (fait) "
              + f"{stoptime-starttime:.3f}s")

    if showtimes:
        print(f"Calcul des blocs d'entiers...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    # Calcul des blocs d'entiers suivant longueur et nombre d'occurrences.
    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level > 2:
        block3 = blocks[2]
    if level > 3:
        block4 = blocks[3]
    maxblock = blocks[-1]
    # ATTENTION!
    # COMPARED TO FEB 2024 VERSION WE SHIFT BY +1 ALL INTEGERS IN
    # SUBLISTS OF maxblock. This is to avoid having to use n+1
    # afterwards for inverse power sums.
    maxblockshifted = [[ n + 1  for n in L] for L in maxblock]
    # NOTA BENE: maxblockshifted will have as last element an empty []
    #            if d=0
    #            This empty [] will not cause problems for the sum()'s
    #            such as sum(1/Rmax(x) for x in maxblockshifted[i])

    # calcul parallèle des beta (sommes d'inverses de puissances).
    # For comments, see irwin().
    # Pay attention though that _lesbetas_par_nb_occurrences
    # produces here beta's which are sums of 1/(n+1)**(m+1)'s for certain
    # n's whereas in irwin() it was sums of 1/n**(m+1).
    # Hence the word "shifted" and usage of maxblockshifted.
    if showtimes:
        stoptime = time.perf_counter()
        print("{:.3f}s".format(stoptime - starttime))

        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
        _lesbetas_par_nb_occurrences = _v5_map_beta_withtimes(Mmax,
                                                              IndexToR,
                                                              maxblockshifted)
    else:
        _lesbetas_par_nb_occurrences = _v5_map_beta_notimes(Mmax,
                                                            IndexToR,
                                                            maxblockshifted)

    # This returns the list L0 such that L0[m] is the sum of the 1/(n+1)**(m+1)
    # where n has level digits and none of them is d.
    # If showtimes is True it prints timings.
    lesbetas_maxblockshifted0 = _lesbetas_par_nb_occurrences(0)
    # idem for 1 occurrence of d
    if k >= 1:
        lesbetas_maxblockshifted1 = _lesbetas_par_nb_occurrences(1)
    # idem for 2 occurrences of d
    if k >= 2:
        lesbetas_maxblockshifted2 = _lesbetas_par_nb_occurrences(2)
    # idem for 3 occurrences of d
    if (k >= 3) and (level > 2):
        lesbetas_maxblockshifted3 = _lesbetas_par_nb_occurrences(3)
    # idem for 4 occurrences of d
    if (k >= 4) and (level > 3):
        lesbetas_maxblockshifted4 = _lesbetas_par_nb_occurrences(4)

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
            S = sum(1/Rmax(x) for x in block1[0])
        elif j == 1:
            if d != 0:
                S = 1/Rmax(d)

        if verbose:
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(S)

        if 2 < level:
            if j <= 2:
                S += sum(1/Rmax(x) for x in block2[j])
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(S)

        if 3 < level:
            if j <= 3:
                S += sum(1/Rmax(x) for x in block3[j])
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(S)

        # ATTENTION
        #
        # So far the S value is the same as for irwin(). But for the
        # positive series, the next contribution in the Burnol
        # series is b times the sum of the 1/(n+1) (not 1/n) where
        # the n's have exactly "level" digits and *at most* j
        # occurrences of digit d.  This is why we have the
        # maxblockshifted[i] here.
        S += b * (sum(sum(1/Rmax(x) for x in maxblockshifted[i])
                      for i in range(1 + min(j,level))))

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(S)

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série positive")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.perf_counter()

        # We now compute the Burnol series which is the positive series
        # as in equation (4) (Theorem 4) of arXiv:2402.09083.
        # We start with the smallest term contributing to the series

        # Each integer n with level digits and no occurrence of d contributes
        # v_{j;m}/(n+1)**(m+1).
        Rm = IndexToR[-1]
        bubu = touslescoeffs[Mmax][j] * lesbetas_maxblockshifted0[Mmax]
        # Each integer n with level digits and 1 occurrence of d contributes
        # v_{j-1;m}/(n+1)**(m+1).
        if j >= 1:
            bubu += touslescoeffs[Mmax][j-1] * lesbetas_maxblockshifted1[Mmax]
        # Each integer n with level digits and 2 occurrences of d contributes
        # v_{j-2;m}/(n+1)**(m+1).
        if j >= 2:
            bubu += touslescoeffs[Mmax][j-2] * lesbetas_maxblockshifted2[Mmax]
        # Each integer n with level digits and 3 occurrences of d contributes
        # v_{j-3;m}/(n+1)**(m+1).
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[Mmax][j-3] * lesbetas_maxblockshifted3[Mmax]
        # Each integer n with level digits and 4 occurrences of d contributes
        # v_{j-4;m}/(n+1)**(m+1).  This can happen here only if level is 4,
        # as we have not implemented higher levels.  And if d is non zero,
        # there is one contribution, if d is zero there are none.
        if (level == 4) and (j >= 4):
            bubu += touslescoeffs[Mmax][j-4] * lesbetas_maxblockshifted4[Mmax]

        if verbose:
            lastterm = bubu  # The Feb 2024 version had a bug here in this
                             # branch printing some info, it had kept -bubu
                             # copied-pasted from irwin() code.
                             # This had no incidence on final result.
            if float(lastterm) == 0.:
                u, E = _v5_shorten_small_real(lastterm)
                print(f"The {Mmax}th term is about {u:f} 10^{E} i.e. ",
                      end = "", flush = True)
            else:
                print(f"The {Mmax}th term is about {lastterm:.3e} i.e. ",
                      end = "", flush = True)

        # COMPUTATION OF THE MAIN SERIES BUILDING UP FROM SMALLEST TERMS
        # Rm will be the RealField. When m decreases Rm changes from time to
        # time regularly and automatically to use more bits.
        for m in range(Mmax-1, 0, -1):  # last one is m=1
            Rm = IndexToR[m]
            # Extend the partial sum obtained to a RealField using
            # higher precision to prepare addition of a new term,
            # which is computed with higher precision.
            bubu = Rm(bubu)
            # and add the new term
            # (if previous step had been  skipped, the value here would
            #  be coerced to the lower precision, hence a completely wrong
            #  final result).
            #
            # See comments above about the contributions 1/(n+1)**(m+1)
            # for integers n having level digits, depending on the count
            # of d's.
            bubu += touslescoeffs[m][j] * lesbetas_maxblockshifted0[m]

            if j >= 1:
                bubu += touslescoeffs[m][j-1] * lesbetas_maxblockshifted1[m]
            if j >= 2:
                bubu += touslescoeffs[m][j-2] * lesbetas_maxblockshifted2[m]
            if (level > 2) and (j >= 3):
                bubu += touslescoeffs[m][j-3] * lesbetas_maxblockshifted3[m]
            if (level == 4) and (j >= 4):
                bubu += touslescoeffs[m][j-4] * lesbetas_maxblockshifted4[m]

        if showtimes:
            stoptime = time.perf_counter()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = S + bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = _v5_shorten_small_real(ratio)
                print("%.3f 10^%s of the total." % (u, E))
            else:
                print("%.3e of the total." % ratio)

            print("La somme de m=1 à %s vaut" % Mmax)
            print(bubu)  # Feb 2024 version printed -bubu by mistake

        if all:
            Sk.append(S)

    if all:
        for j in range(k+1):
            print(f"(k={j}) {Rfinal(Sk[j])}")

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    return Rfinal(S)


if __name__ == "__main__":
    print(f"""
Hello, this file {__filename__} provides two functions irwin()
and irwinpos().  Use help(irwin) or help(irwinpos) for help.

This is version {__version__} of {__date__}.

General information is also available in the irwin_v5_docstring
variable.

The variable "maxworkers" sets ncpus for @parallel usage. You
can set it and then reload.
{maxworkersinfostring}
"""
          )
