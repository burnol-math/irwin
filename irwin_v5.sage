# -*- mode: python ; coding: utf-8; -*-

# irwin_v5.sage
# Use via load("irwin_v5.sage") in sage interactive mode

__version__  = "1.5.7"
__date__     = "2025/05/16"
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
    parallel with m=M+1, M+2, ...., M+n with a=1, 2, ..., q.
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
    return list(_v5_beta_aux(m, IR[m], nblock)
                for m in range(start, end, maxworkers))


def _v5_map_beta_notimes(Mmax, IndexToR, maxblock):
    """Auxiliary for sharing code between irwin() and irwinpos().

    TODO: explain what it does.
    """
    extra = maxworkers - ( Mmax % maxworkers )
    if extra < maxworkers:
        def map__v5_beta(j):
            """TODO: add some docstring.
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
            """TODO: add some docstring.
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
    """Auxiliary for sharing code between irwin() and irwinpos().
    TODO: explain what it does.
    """
    # We want to display some visual sign of progress.
    # Find the largest multiple of maxworkers at most 1000,
    # do something reasonable if maxworkers is big
    q = max(1000 // maxworkers, 32)
    mSize = q * maxworkers
    def map__v5_beta(j):
        """TODO: add some docstring.
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
    """Organize integers according to nb of digits and d count
    """
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
        print("Calcul des blocs initiaux et des gammas...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level >2:
        block3 = blocks[2]
        if level >3:
            block4 = blocks[3]
    maxblock = blocks[-1]
    # NOTA BENE: maxblock will have as last element an empty [] if d=0
    #            This empty [] will not cause problems for the sum()'s
    #            such as sum(1/Rmax(x) for x in maxblock[i])

    # lesgammas[j] is only ever needed to compute a u_{k;m} for m
    # at least equal to j. It is used only in a sum with
    # non-negative contributions, there are no subtractions. So we
    # only need it to the precision needed for u_{k;m}
    # itself. There is however a cumulative error effect in a sum
    # with m contributions, but it should be absorbed in the
    # safety cushion we have in place in evaluating the needed
    # precision.  Perhaps if k is very large, this safety cushion
    # could prove defective.
    lesgammas = [ bmoinsun ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lesgammas.append(Rj(sum(a**j for a in block1[0])))

    # Those are only needed for k>O.  Same remark as for
    # lesgammas[j] relative to the precision to use.
    if k > 0:
        lespuissancesded = [ 1 ]
        for j in range(1, Mmax+1):
            Rj = IndexToR[j]
            lespuissancesded.append(Rj(d**j))
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
    Q, R = divmod(Mmax - 1, maxworkers)
    for P in range(Q):
        useparallel = _v5_para_recurrence(m, maxworkers, useparallel)
        m += maxworkers
    # Ici on va invoquer une procédure parallélisée avec <maxworkers.
    if R > 0:
        _ = _v5_para_recurrence(m, R, useparallel)

    if showtimes:
        stoptime = time.perf_counter()
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} "
              + f"Fini! En tout : {stoptime-starttime:.3f}s")

    # Calcul parallèle des beta (sommes d'inverses de puissances).
    if showtimes:
        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
        map__v5_beta = _v5_map_beta_withtimes(Mmax, IndexToR, maxblock)
    else:
        map__v5_beta = _v5_map_beta_notimes(Mmax, IndexToR, maxblock)

    lesbetas_maxblock0 = map__v5_beta(0)

    if k >= 1:
        lesbetas_maxblock1 = map__v5_beta(1)

    if k >= 2:
        lesbetas_maxblock2 = map__v5_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblock3 = map__v5_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblock4 = map__v5_beta(4)

    # Boucle qui évalue également la série pour les j<k (si all = True).
    Sk = []

    for j in range(0 if all else k, k+1):
        if showtimes:
            print(f"Calcul de l'approximation principale avec k={j}...",
                  end = ' ', flush = True)
            starttime = time.perf_counter()

        # Calcul de la série alternée de Burnol.

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

        # WE START WITH THE SMALLEST TERM CONTRIBUTING TO THE SERIES
        # Its sign will be set later.

        Rm = IndexToR[-1]

        # Compared to 2024 version touslescoeffs has its two indices permuted
        bubu = touslescoeffs[Mmax][j] * lesbetas_maxblock0[Mmax]

        if j >= 1:
            bubu += touslescoeffs[Mmax][j-1] * lesbetas_maxblock1[Mmax]
        if j >= 2:
            bubu += touslescoeffs[Mmax][j-2] * lesbetas_maxblock2[Mmax]
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[Mmax][j-3] * lesbetas_maxblock3[Mmax]
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
            bubu = Rm(-bubu)
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
        print("Calcul des blocs initiaux et des gammas ...",
              end = ' ', flush = True)
        starttime = time.perf_counter()

    blocks = _v5_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level >2:
        block3 = blocks[2]
        if level >3:
            block4 = blocks[3]
    maxblock = blocks[-1]
    # COMPARED TO FEB 2024 VERSION WE SHIFT BY +1 ALL INTEGERS IN
    # SUBLISTS OF maxblock. This is to avoid having to use n+1
    # afterwards for inverse power sums.
    maxblockshifted = [[ n + 1  for n in L] for L in maxblock]
    # NOTA BENE: maxblockshifted will have as last element an empty []
    #            if d=0
    #            This empty [] will not cause problems for the sum()'s
    #            such as sum(1/Rmax(x) for x in maxblockshifted[i])


    # ATTENTION que la série positive a des récurrences avec b-1-d
    # à la place de d
    Aprime = [i for i in range(b)]
    dprime = b - 1 - d
    Aprime.remove(dprime)

    # See comments in irwin() for the precision used.
    lesgammasprime = [ bmoinsun ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lesgammasprime.append(Rj(sum(a**j for a in Aprime if a != 0)))

    # Same reasoning for the precision.
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
    # Attention to this b * b  extra needed for the  v_{0;1}.
    c1 = [ (b * b + lesgammasprime[1] * Rmax(b)) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(( (lesgammasprime[1] + dprime) * Rmax(b) + c1[-1])/Rmax(b * b - bmoinsun))
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

    # calcul parallèle des beta (sommes d'inverses de puissances).
    if showtimes:
        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
        map__v5_beta = _v5_map_beta_withtimes(Mmax, IndexToR, maxblockshifted)
    else:
        map__v5_beta = _v5_map_beta_notimes(Mmax, IndexToR, maxblockshifted)

    lesbetas_maxblockshifted0 = map__v5_beta(0)

    if k >= 1:
        lesbetas_maxblockshifted1 = map__v5_beta(1)

    if k >= 2:
        lesbetas_maxblockshifted2 = map__v5_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblockshifted3 = map__v5_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblockshifted4 = map__v5_beta(4)

    # Boucle qui évalue également la série pour les j<k (si all = True).
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

        # Attention à emploi de maxblockshifted ici.
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

        # WE START WITH THE SMALLEST TERM CONTRIBUTING TO THE SERIES
        # Its sign will be set later.

        # Attention à emploi de maxblockshifted ici.
        Rm = IndexToR[-1]
        bubu = touslescoeffs[Mmax][j] * lesbetas_maxblockshifted0[Mmax]

        if j >= 1:
            bubu += touslescoeffs[Mmax][j-1] * lesbetas_maxblockshifted1[Mmax]
        if j >= 2:
            bubu += touslescoeffs[Mmax][j-2] * lesbetas_maxblockshifted2[Mmax]
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[Mmax][j-3] * lesbetas_maxblockshifted3[Mmax]
        if (level == 4) and (j >= 4):
            bubu += touslescoeffs[Mmax][j-4] * lesbetas_maxblockshifted4[Mmax]

        if verbose:
            lastterm = bubu  # the Feb 2024 version had a bug here
                             # doing -bubu if Mmax was odd
                             # (as copy pasted from alternating series code)
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

        # Attention à emploi de maxblockshifted ici.
        for m in range(Mmax-1, 0, -1):  # last one is m=1
            Rm = IndexToR[m]
            # extend partial sum to higher precision
            bubu = Rm(bubu)
            # and add new term (if previous step is skipped, the result will get
            #                   coerced to lower precision, hence ultimately the
            #                   final result is completely wrong).
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
