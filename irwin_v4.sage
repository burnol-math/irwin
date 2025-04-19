# -*- mode: python ; coding: utf-8; -*-

# irwin_v4.sage
# Use via load("irwin_v4.sage") in sage interactive mode

__version__  = "1.4.7"
__date__     = "2025/04/18"
__filename__ = "irwin_v4.sage"

irwin_v4_docstring = """
This file is an evolution of irwin.sage as available at

    https://arxiv.org/src/2402.09083v5/anc

with three main enhancements relevant to targeting 1000 or
more decimal digits of precision:

- Use of decreasing precision for higher terms in the series.
  This is done with a granularity of PrecStep bits.  It is an
  optional parameter to irwin() and irwinpos() which defaults to
  500.

- There is no pre-computation of the Pascal triangle up to N=1000
  anymore.  Only two rows of the Pascal triangle are kept at any
  given time in memory.  This reduces the memory footprint by a
  factor of N/6: the number of bytes needed goes down from about
  N^3/(96 log(2)) (this takes into account a factor 2 from using
  the symmetry binom(n,j)=binom(n,n-j)) to only N^2/(16 log(2))
  (1byte=8bits).

  For N=10000, N^3/(96 log(2)) is about 1.5e10 i.e. 14GiB or 15GB.
  But N^2/(16 log(2)) is only about 9e6 i.e. 8.6MiB or 9MB.  For
  higher N's the N/6 factor makes the difference between the code
  being usable or not usable...

- Parallelization (using @parallel(ncpus=8) per default) of
  some parts of the algorithm, mainly the computation of the
  beta coefficients.

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
def docstring_parameter(*sub):
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj
    return dec

irwin_fn_docstring = """\
    :param int b: the integer base
    :param int d: the digit
    :param int k: the number of occurrences
    :param int nbdigits: (optional, default 34)
        The whished-for number of decimal digits for the result.
    :param int level: (optional, default 3)
        The level must be 2, 3 or 4.
    :param int PrecStep: (optional, default ``500``)
        Terms of the series are computed with a RealField
        of evolving precision, which differs from the maximal
        precision by a suitable multiple of PrecStep.
    :param bool all: (optional, default ``False``)
        If ``True``  print all irwin sums for ``j`` occurrences
        with ``j`` from ``0`` to ``k``.
    :param bool showtimes: (optional, default ``False``)
        Whether to print out timings for various steps.
    :param bool verbose: (optional, default ``False``)
        Whether to print the values of intermediate contributions
        to the final value, in particular to confirm enough
        terms of the series were used.
    :param int Mmax: (optional, default ``-1``)
        If not ``-1`` the number of terms of the series to use.
        Use only if the auto-choice is excessive (this will be
        the case for k=0 and d=1, and to a lesser extent when
        d=b-1).  Use ``verbose=True`` to check how many terms
        are used by default.

    :rtype: :class:`sage.rings.real_mpfr.RealNumber`
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b

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

    sage: irwin(10, 9, 4, 52, all=True)
    (k=0) 22.92067661926415034816365709437593191494476243699848
    (k=1) 23.04428708074784831967594930973617482538959203064774
    (k=2) 23.02604026596124378845022249787272342108112267542086
    (k=3) 23.02585299837244431714290384468012275518705238435290
    (k=4) 23.02585095265829261377053973815542996035002267989413
    23.02585095265829261377053973815542996035002267989413
"""

import time

# print(f"Sage voit {sage.parallel.ncpus.ncpus()} CPUs, cela est-il normal ?")
# voir https://ask.sagemath.org/question/42439/parallel-how-to-use-all-cpus/
# may je dois expérimenter avec SAGE_NUM_THREADS si ça existe encore.
try:
    maxworkers_test = maxworkers
    maxworkersinfostring = ("maxworkers variable already existed "
                            f"and has value {maxworkers}.")
except NameError:
    maxworkers = 8
    maxworkersinfostring = ("maxworkers variable has been created "
                            "and assigned value 8.")

assert maxworkers < 1000, f"Sorry maxworkers={maxworkers} must be less than 1000"

@parallel(ncpus=maxworkers)
def para_um(a, P, G, T, R, j, m):
    return sum(P[i]*R(G[i])*R(T[m-i][j]) for i in range(a,m+1,maxworkers))

def para_beta_aux(m, R, nblock):
    return sum(1/R(n ** (m+1)) for n in nblock)

# legacy comment, kept in case
#   Si p_iter="multiprocessing"
#   ValueError: Cannot pickle code objects from closures
# TODO: check efficiency of way of passing the RealField.
#       earlier version passed directly a RealField
#       here we pass the whole IndexToR, but as
#       this is called now only exactly maxworkers
#       times, it should be rather efficient. I don't
#       have time now, and will not have access to
#       multi-core computer for a while.
@parallel(ncpus=maxworkers)
def para_beta(mlist, IR, nblock):
    return list(para_beta_aux(m, IR[m], nblock) for m in mlist)

nbguardbits = 12

def _shorten_small_real(rr):
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


def _setup_realfields(nbdigits, PrecStep, b, level, Mmax=-1):
    """Preparation of an array mapping index to a RealField

    NOTE BENE: with d=b-1, there is an additional factor
               ((b-2)/(b-1))**m in u_{0;m} which is not taken
               into account here, so we end up using more terms
               than needed in this case.  For example at level =
               2, computing the classic Kempner sum for 500
               digits will use about 500 terms but 475 would
               have been enough as (8/9)**475 is about 5e-25.
               For k>0 and increasing, the effect diminishes.
               (One can use verbose=True to examine the size of
               the smallest kept term).

               For the positive series, there is no such factor
               but the beta(m+1) is smaller for level=2 roughly
               by (b/(b+1))**m (for d=1 this will be rather with
               b/(b+1/2), and d=0 also has another ratio) so we
               have a similar phenomenon of using too many terms
               but for other reasons.  For level=3 (which is the
               default) effect will be lower than for level=2.

               With k=0 (only) and excluded digit 1, we use way
               too many terms due to roughly a 2**-m factor not
               taken into account.  User can set Mmax
               explicitly to their own choosing.

               Here we choose Mmax according to an analysis
               which is supposed to work for all k's and d's
               and level>1.
    """

    # Try to guarantee we will have nbdigits radix 10 digits correct.
    nbbits_final = ceil((nbdigits+1)*log(10,2))
    Rfinal = RealField(nbbits_final)
    # I experimented a bit with using only multiples of 32bits.
    # But for say about digits, it causes the risk of doing
    # computations corresponding to perhaps say 9 more decimal
    # digits in extra of the already about 4 from the guard bits,
    # inducing use of more terms of the series.  Perhaps there is
    # an advantage for tens of thousands of digits when
    # computation times may run into hours, but it is lengthy to
    # test out...
    #
    # The ceiling of maximal precision to a multiple of 32bits was
    # done if PrecStep was multiple of 8. Now dropped.
    # This was mainly done at a time where I experimented with
    # a PrecStep less than 50, but finally it looks preferable
    # to set it to a higher value.
    if True:
        nbbits = nbbits_final + nbguardbits
    else:
        nbbits = 32 * ceil((nbbits_final + nbguardbits)/32)
    R = RealField(nbbits)

    # We estimate how many terms are needed according to level and
    # targeted final precision. This estimate relies on two things:
    # - a priori lower bound for the final result
    # - a priori upper bounds for the terms of the series.

    # The u_{k;m}'s are bounded above by b. In fact they are
    # bounded above by b/(m+1) but we don't have this upper bound
    # for the v_{k;m}'s, on the other hand the beta's are smaller
    # for them as they use inverse powers of n+1 not n.

    # We always take level at least 2.  The crucial upper bound
    # comes from the beta(m+1)'s.  We can estimate that beta(m+1)
    # is bounded above by (1/b**(l-1) + 1/m)/b**((l-1)m) with here
    # l=level.  For k=0 and d=1, we can improve that to
    # (1/2b**(l-1)+1/m)/2**m/b**((l-1)m).

    # We know from Farhi's theorem that blog(b) is a lower bound
    # of total sum if k>0 or if k=0 and d=0.  If k=0 and d>0 we
    # have Proposition 6 in my Irwin paper which says S > b*log(b)
    # -b*log(1+1/d). Worst case is with d=1, giving b*log(b/2). For
    # d at least 2 the lower bound improves to b*log(2b/3). Using
    # only that the u_{k;m}'s and v_{k;m}'s are bounded above by
    # b, we will need to compare log(b/2) with
    # (1/2b**(l-1)+1/m)/2**m, and log(2b/3) with 1/b**(l-1) + 1/m.
    # We always take l at least 2.  So it is log(b/2) versus
    # (1/(2b)+1/m)/2**m or log(2b/3) versus 1/b + 1/m.
    #
    # Testing the worst case m=1, we have log(b/2)>(1/(2b)+1)/2
    # starting at b=4, and log(2b/3)>1/b + 1 for b at least 5.
    #
    # For m=2, we have log(b/2)>(1/(2b)+1/2)/4 for b at least 3
    # and log(2b/3)>1/b+1/2 for b at least 4.  But for b=3, the
    # difficulty with m=2 is for for k=0 and d neither 0 (as
    # log(3)>1>1/3+1/2) or 1, hence d=2. One finds numerically
    # S>2.682, and S/3>0.894>1/3+1/2.
    #
    # Only remains to check for b=2.  The only Irwin sum not
    # >2log(2) is with k=0 and d=1, it is the empty sum with value
    # zero. Use irwin(2,1,0) or irwinpos(2,1,0) at your own
    # risk...  (it seems to work fine).
    #
    # For level=2, beta(m+1) for b=2 is bounded above by
    # 1/2**(m+1)+1/3**(m+1) and u_{k;m} by 2/(m+1).  We want to
    # compare this with 2 log(2)/2**m. It is less already for m=1.
    #
    # For the positive series we have 2(1/3**(m+1)+1/4**(m+1)) to
    # compare with 2log(2)/2**m. Again it is less already for m=1.
    #
    # For level>2 and the alternating series a more precise upper
    # bound of the m th term is (1/b**(l-1)+1/m)b/(m+1) times
    # b**(-(l-1)m) so we compare (1/4+1/m)/(m+1) with log(2).  And
    # already for m=1 it is smaller. Still for b=2 and l>2, for
    # the positive series we can bound above the mth term by
    # (1/(b**(l-1)+1)+1/m) * (b**(l-1)/(b**(l-1)+1))**m * b times
    # b**(-(l-1)m) using only v_{k;m}<= b. And we need to compare
    # with 2 log(2). So we check if (1/5+1/m)*(4/5)**m is less
    # than log(2). This is true for m at least 2.

    # In conclusion it is always true that for m at least 2 the m
    # th term of the series contributes a fraction less than
    # b**(-(l-1)m) of the total sum (we could make a detailed
    # examination for m=1, but drop it). For the alternative
    # series this gives a bound for the total error made by
    # neglecting all terms starting with the m th one; for the
    # positive series, we have to take into account an extra
    # factor of b**(l-1)/(b**(l-1)-1) which is at most 2.
    # We don't worry about this 2.

    # We only need to take into account the m th term of the series
    # with a precision equal to nbbits - (l-1)*m*log(b,2). (the first
    # one will be computed with full precision).
    # We choose the number of terms Mmax such that it is the last
    # with nbbits - (l-1)*Mmax*log(b,2) >= nbguardbits/2
    _Mmax = floor((nbbits - nbguardbits/2)/(level-1)/log(b,2))

    # The mth term only needs to be computed with precision
    # nbbits - (l-1)*m*log(b,2) but we will use precisions
    # of the type nbbits - j T with T = PrecStep.
    # So we have nbbits -jT >= nbguardbits/2
    NbOfPrec = 1 + floor((nbbits - nbguardbits/2)/PrecStep)

    # We add an overhead or pre-creating all used RealFields.
    # Maybe a bit stupid.
    LesReels = [R]
    for j in range(1, NbOfPrec):
        LesReels.append(RealField(nbbits - j * PrecStep))

    # We need to map m to j in such a way that j is the largest such
    # that nbbits - jT >= nbbits - (l-1)*m*log(b,2)
    # So j is floor((l-1)*log(b,2)*m/T)
    #
    # We add the overhead of a list mapping m to the RealField
    # it will use.  This is to avoid having to do dynamically
    # j = floor((l-1)*log(b,2)*m/T) and this all is probably
    # very silly.
    #
    # To do this we need to know which m will use a given j.
    # They verify j<= (l-1)*log(b,2)*m/T < j+1
    # j*T/(l-1)/log(b,2)<= m < (j+1)*T/(l-1)/log(b,2)
    # ceil(j*T/(l-1)/log(b,2))<= m < ceil((j+1)*T/(l-1)/log(b,2)).
    #
    # If T/(l-1)/log(b,2)<1 it is not guaranteed that a given j will
    # have associated m's. But this is not a problem.
    IndexToPrecRatio = PrecStep / (level - 1) / log(b,2)  # not an integer !
    IndexToR = []
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


def _setup_blocks(b, d, level):
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


def _comp_um_print_timeinfo(single, multi, para, wrkrs, m):
    """Auxiliary shared between irwin() and irwinpos()
    """
    if multi < single:
        if para:
            print(f"\nPoursuite car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle (m={m})", end="")
        else:
            print(f"\nBasculement car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle (maxworkers={wrkrs}, m={m})", end="")
    else:
        if para:
            print(f"\nOn quitte car {single:.3f}s<{multi:.3f}s"
                  f" l'exécution parallèle (m={m}) ", end="")
        else:
            print(f"\nPas utile ({single:.3f}s<{multi:.3f}s)"
                  f" d'exécuter en parallèle (m={m}) ", end="")


@docstring_parameter(irwin_fn_docstring)
def irwin(b, d, k,
          nbdigits=34,
          level=3,
          PrecStep=500,
          all=False,
          showtimes=False,
          verbose=False,
          Mmax=-1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise l'algorithme de Burnol, série alternée de niveau 2, 3 ou 4.

    {0}
    """

    assert 1 < level <= 4, "Le niveau (level) doit être 2 ou 3 ou 4"

    # SURTOUT NE PAS FAIRE if type(b) == type(1) !!!!
    # ça marche en pour des inputs directs mais pas pour "for b in range(B)"
    # J'avais bêtement ajouté cette erreur le samedi 24 février 2024 v2 sur arXiv
    assert b > 1, "%s doit être au moins 2" % b
    bmoinsun = b - 1

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    if showtimes:
        print("Préparation des RealField...", end=" ", flush=True)
        starttime = time.time()

    (nbbits,
     nbbits_final,
     R,
     Rfinal,
     Mmax,
     IndexToR,
     NbOfPrec) = _setup_realfields(nbdigits, PrecStep, b, level, Mmax)

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} RealField(s) de précision maximale {nbbits},")
        print(f"décrémentée par multiples de {PrecStep}")

    if showtimes:
        print("Calcul des blocs initiaux et des gammas...",
              end = ' ', flush = True)
        starttime = time.time()

    blocks = _setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level >2:
        block3 = blocks[2]
        if level >3:
            block4 = blocks[3]
    maxblock = blocks[-1]
    # NOTA BENE: maxblock will have as last element an empty [] if d=0
    #            This empty [] will not cause problems for the sum()'s
    #            such as sum(1/R(x) for x in maxblock[i])

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
    lespuissancesded = [ 1 ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lespuissancesded.append(Rj(d**j))

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))

    if showtimes:
        print(f"Calcul des coefficients u_{{j;m}} pour 0<=j<={k}...",
              end = ' ', flush = True)
        starttime = time.time()

    # calcul récursif des moments
    # Je pense (2025) que j'avais conclu qu'il fallait mieux
    # calculer les puissances comme des entiers exactement.
    # In order to have to evaluate each Pascal triangle row only once, we
    # use a slightly modified syntax compared to the 2024 version, instead
    # of having touslescoeffs = [ [ the u_{0,m}'s ], [ the u_{1,m}'s ], ... ]
    # it will be touslescoeffs = [[u_{0,0}, u_{1,0}, ..., u_{k,0}],
    #                             [u_{0,1}, u_{1,1}, ..., u_{k,1}], ... ]
    touslescoeffs = [ [R(b)] * (k+1) ]
    c1 = [ lesgammas[1] * R(b) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(((lesgammas[1] + d) * R(b) + c1[-1])/R(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRow = [1, 1]
    useparallel = False
    if Mmax > 511 and showtimes and verbose:
        print("\n    L'algorithme teste tous les 512 coefficients s'il est bénéfique")
        print("    d'utiliser la procédure décorée par @parallel.  Ce test n'est fait")
        print("    que pour le calcul des u_{0;m}.  Le premier nombre indiqué sera")
        print("    la durée du calcul normal, le second avec @parallel, pour u_{0;m}.",
              end = "", flush = True)

    for m in range(2, Mmax+1):
        prevPascalRow = PascalRow.copy()
        PascalRow = [1]
        halfm = m//2
        PascalRow.extend([prevPascalRow[j-1]+prevPascalRow[j] for j in range(1, halfm)])
        if not (m&1):
            L = PascalRow.copy()
        PascalRow.append(prevPascalRow[halfm-1]+prevPascalRow[halfm])
        if m&1:
            PascalRow.extend(reversed(PascalRow))
        else:
            PascalRow.extend(reversed(L))

        Rm = IndexToR[m]
        if m&511:
            if useparallel:
                results = list(para_um(((a, PascalRow, lesgammas, touslescoeffs, Rm, 0, m,)
                                       for a in range(1, maxworkers + 1))))
                partial_sums = [result[1] for result in results]
                # De toute façon on ne peut pas faire que ce soit exactement
                # le même calcul qu'en non-parallèle puisque les termes ont
                # déjà été additionnés modulo maxworkers.  Faire sorted()
                # ici n'est pas nécéssaire, et ne semble pas avoir d'intérêt
                # particulier d'ailleurs.  On pourrait faire sorted(partial_sums)
                # mais bon.
                cm = [ sum(partial_sums) / Rm(b**(m+1) - bmoinsun) ]
            else:
                cm = [ sum(PascalRow[i]*Rm(lesgammas[i])*Rm(touslescoeffs[m-i][0])
                           for i in range(1,m+1)) / Rm(b**(m+1) - bmoinsun) ]
        else:
            localstarttime = time.time()
            results = list(para_um(((a, PascalRow, lesgammas, touslescoeffs, Rm, 0, m,)
                                    for a in range(1, maxworkers + 1))))
            partial_sums = [result[1] for result in results]
            x = [ sum(partial_sums) / Rm(b**(m+1) - bmoinsun) ]
            multitime = time.time() - localstarttime

            localstarttime = time.time()
            cm = [ sum(PascalRow[i]*Rm(lesgammas[i])*Rm(touslescoeffs[m-i][0])
                       for i in range(1,m+1)) / Rm(b**(m+1) - bmoinsun) ]

            singletime = time.time() - localstarttime

            if showtimes:
                _comp_um_print_timeinfo(singletime, multitime,
                                        useparallel, maxworkers, m)

            useparallel = (multitime < singletime)

        if useparallel:
            for j in range(1,k+1):
                results = list(para_um(((a, PascalRow, lesgammas,
                                            touslescoeffs, Rm, j, m,)
                                        for a in range(1, maxworkers + 1))))
                x = sum(p[1] for p in results)
                x += cm[-1]
                results = list(para_um(((a, PascalRow, lespuissancesded,
                                            touslescoeffs, Rm, j-1, m,)
                                        for a in range(1, maxworkers + 1))))
                x += sum(p[1] for p in results)
                x = x / Rm(b**(m+1) - bmoinsun)
                cm.append(x)
        else:
            for j in range(1,k+1):
                cm.append((sum(PascalRow[i]*Rm(lesgammas[i])*Rm(touslescoeffs[m-i][j])
                               for i in range(1,m+1))
                           + cm[-1]
                           + sum(PascalRow[i]*Rm(lespuissancesded[i])*
                                 Rm(touslescoeffs[m-i][j-1])
                                 for i in range(1,m+1))
                           )
                          / Rm(b**(m+1) - bmoinsun)
                          )

        touslescoeffs.append(cm.copy())

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))

    # calcul parallèle des beta (sommes d'inverses de puissances)
    if showtimes:
        print(f"Calcul parallélisé des beta(m+1) avec maxworkers={maxworkers}")
        # We want to display some visual sign or progress
        # Find the largest multiple of maxworkers at most 1000
        # HOWEVER THIS COSTS IN CODE SIMPLICITY AND EFFICIENCY
        # ATTENTION maxworkers must be < 1000...
        mSize = (1000 // maxworkers) * maxworkers
        q, r = divmod(mSize, maxworkers)
        I = 0
        indices = []
        for i in range(r):
            indices.append(I)
            I += q + 1
        for i in range(maxworkers - r):
            indices.append(I)
            I += q
        # assert I == mSize, "check your math"
        indices.append(I)
        def map_para_beta(j):
            print(f"({j} occ.)...", end = ' ', flush = True)
            starttime = time.time()
            mend = 1
            L = [0]
            for i in range(Mmax//mSize):
                mbegin = 1 + i*mSize
                mend   = mbegin + mSize
                # TODO: find a way to not "instantiate" m range and keep
                #       it as a generator but I have to become more
                #       knowledgeable in Python and impact will be minor
                #       anyhow
                mrange = list(range(mbegin, mend))
                inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblock[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(para_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{mend}", end = " ", flush= True)

            if mend < Mmax+1:
                mrange = list(range(mend, Mmax + 1))
                qlast, rlast = divmod(len(mrange), maxworkers)
                I = 0
                indiceslast = []
                for i in range(rlast):
                    indiceslast.append(I)
                    I += qlast + 1
                for i in range(maxworkers - rlast):
                    indiceslast.append(I)
                    I += qlast
                indiceslast.append(I)
                inputblocks = [(mrange[indiceslast[i]:indiceslast[i+1]], IndexToR,
                                maxblock[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(para_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{Mmax+1} (done)", end = " ", flush=True)
            stoptime = time.time()
            print("{:.3f}s".format(stoptime - starttime))
            return L
    else:
        q, r = divmod(Mmax, maxworkers)
        I = 0
        indices = []
        for i in range(r):
            indices.append(I)
            I += q + 1
        for i in range(maxworkers - r):
            indices.append(I)
            I += q
        indices.append(I)
        def map_para_beta(j):
            L = [0]
            mrange = list(range(1, Mmax + 1))
            inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblock[j])
                           for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(para_beta(inputblocks)))]
            L.extend(sum([result for result in results_1], []))
            return L

    lesbetas_maxblock0 = map_para_beta(0)

    if k >= 1:
        lesbetas_maxblock1 = map_para_beta(1)

    if k >= 2:
        lesbetas_maxblock2 = map_para_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblock3 = map_para_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblock4 = map_para_beta(4)

    # boucle pour évaluer également les j < k si all = True
    Sk = []

    for j in range(0 if all else k, k+1):
        if showtimes:
            print(f"Calcul de l'approximation principale avec k={j}...",
                  end = ' ', flush = True)
            starttime = time.time()

        # calcul de la série alternée de Burnol

        S = 0

        if j == 0:
            S = sum(1/R(x) for x in block1[0])
        elif j == 1:
            if d != 0:
                S = 1/R(d)

        if verbose:
            print("Somme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(S)

        if 2 < level:
            if j <= 2:
                S += sum(1/R(x) for x in block2[j])
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(S)

        if 3 < level:
            if j <= 3:
                S += sum(1/R(x) for x in block3[j])
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(S)

        # jusqu'à j répétitions ; si j >= level, on s'arrête à
        # level répétitions max
        S += b * (sum(sum(1/R(x) for x in maxblock[i])
                      for i in range(1 + min(j,level))))

        if showtimes:
            stoptime = time.time()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(S)

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série alternée")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.time()

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
                u, E = shorten_small_real(lastterm)
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
            stoptime = time.time()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = S - bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = _shorten_small_real(ratio)
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


@docstring_parameter(irwin_fn_docstring)
def irwinpos(b, d, k,
             nbdigits=34,
             level=3,
             PrecStep=500,
             all=False,
             showtimes=False,
             verbose=False,
             Mmax=-1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise algorithme de Burnol, série positive de niveau 2, 3 ou 4.

    {0}
    """

    assert 1 < level <= 4, "Le niveau (level) doit être 2 ou 3 ou 4"

    # SURTOUT NE PAS FAIRE if type(b) == type(1) !!!!
    # ça marche en pour des inputs directs mais pas pour "for b in range(B)"
    # J'avais bêtement ajouté cette erreur le samedi 24 février v2 sur arXiv
    assert b > 1, "%s doit être au moins 2" % b
    bmoinsun = b - 1

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    if showtimes:
        print("Préparation des RealField...", end=" ", flush=True)
        starttime = time.time()

    (nbbits,
     nbbits_final,
     R,
     Rfinal,
     Mmax,
     IndexToR,
     NbOfPrec) = _setup_realfields(nbdigits, PrecStep, b, level, Mmax)

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))
    if verbose:
        print(f"{NbOfPrec} RealField(s) de précision maximale {nbbits},")
        print(f"décrémentée par multiples de {PrecStep}")

    if showtimes:
        print("Calcul des blocs initiaux et des gammas...",
              end = ' ', flush = True)
        starttime = time.time()

    blocks = _setup_blocks(b, d, level)
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
    #            such as sum(1/R(x) for x in maxblockshifted[i])


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
    lespuissancesdedprime = [ 1 ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lespuissancesdedprime.append(Rj(dprime**j))

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))

    if showtimes:
        print(f"Calcul des coefficients v_{{j;m}} pour 0<=j<={k}...",
              end = ' ', flush = True)
        starttime = time.time()

    # calcul récursif des moments
    # Je pense (2025) que j'avais conclu qu'il fallait mieux calculer les
    # puissances comme des entiers exactement.
    # In order to have to evaluate each Pascal triangle row only once, we
    # use a slightly modified syntax compared to the 2024 version, instead
    # of having touslescoeffs = [ [ the v_{0,m}'s ], [ the v_{1,m}'s ], ... ]
    # it will be touslescoeffs = [[v_{0,0}, v_{1,0}, ..., v_{k,0}],
    #                             [v_{0,1}, v_{1,1}, ..., v_{k,1}], ... ]
    touslescoeffs = [ [R(b)] * (k+1) ]
    # attention to b * b  extra in first one
    c1 = [ (b * b + lesgammasprime[1] * R(b)) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(( (lesgammasprime[1] + dprime) * R(b) + c1[-1])/R(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRow = [1, 1]
    useparallel = False
    if Mmax > 511 and showtimes and verbose:
        print("\n    L'algorithme teste tous les 512 coefficients s'il est bénéfique")
        print("    d'utiliser la procédure décorée par @parallel.  Ce test n'est fait")
        print("    que pour le calcul des v_{0;m}.  Le premier nombre indiqué sera")
        print("    la durée du calcul normal, le second avec @parallel, pour v_{0;m}.",
              end = "", flush = True)

    for m in range(2, Mmax+1):
        prevPascalRow = PascalRow.copy()
        PascalRow = [1]
        halfm = m//2
        PascalRow.extend([prevPascalRow[j-1]+prevPascalRow[j] for j in range(1, halfm)])
        if not (m&1):
            L = PascalRow.copy()
        PascalRow.append(prevPascalRow[halfm-1]+prevPascalRow[halfm])
        if m&1:
            PascalRow.extend(reversed(PascalRow))
        else:
            PascalRow.extend(reversed(L))

        Rm = IndexToR[m]

        if m&511:
            if useparallel:
                results = list(para_um(((a, PascalRow, lesgammasprime,
                                         touslescoeffs, Rm, 0, m,)
                                        for a in range(1, maxworkers + 1))))
                partial_sums = [result[1] for result in results]
                # De toute façon on ne peut pas faire que ce soit exactement
                # le même calcul qu'en non-parallèle puisque les termes ont
                # déjà été additionnés modulo maxworkers.  Faire sorted()
                # ici n'est pas nécéssaire, et ne semble pas avoir d'intérêt
                # particulier d'ailleurs.  On pourrait faire sorted(partial_sums)
                # mais bon.
                #  attention to this extra term specific to v_m recurrence
                cm = [ (Rm(b ** (m+1)) + sum(partial_sums)) / Rm(b**(m+1) - bmoinsun) ]
            else:
                #  attention to this extra term specific to v_m recurrence
                cm = [ (Rm(b ** (m+1))
                        + sum(PascalRow[i]*Rm(lesgammasprime[i]) *
                              Rm(touslescoeffs[m-i][0])
                              for i in range(1,m+1))
                        ) / Rm(b**(m+1) - bmoinsun) ]
        else:
            localstarttime = time.time()
            results = list(para_um(((a, PascalRow, lesgammasprime,
                                     touslescoeffs, Rm, 0, m,)
                                    for a in range(1, maxworkers + 1))))
            partial_sums = [result[1] for result in results]
            x = [ (Rm(b ** (m+1)) + sum(partial_sums)) / Rm(b**(m+1) - bmoinsun) ]
            multitime = time.time() - localstarttime

            localstarttime = time.time()
            #  attention to this extra term specific to v_m recurrence
            cm = [ (Rm(b ** (m+1))
                    + sum(PascalRow[i]*Rm(lesgammasprime[i]) *
                          Rm(touslescoeffs[m-i][0])
                          for i in range(1,m+1))
                    ) / Rm(b**(m+1) - bmoinsun) ]

            singletime = time.time() - localstarttime

            if showtimes:
                _comp_um_print_timeinfo(singletime, multitime,
                                        useparallel, maxworkers, m)

            useparallel = (multitime < singletime)

        if useparallel:
            for j in range(1,k+1):
                results = list(para_um(((a, PascalRow, lesgammasprime,
                                            touslescoeffs, Rm, j, m,)
                                        for a in range(1, maxworkers + 1))))
                x = sum(p[1] for p in results)
                x += cm[-1]
                results = list(para_um(((a, PascalRow, lespuissancesdedprime,
                                            touslescoeffs, Rm, j-1, m,)
                                        for a in range(1, maxworkers + 1))))
                x += sum(p[1] for p in results)
                x = x / Rm(b**(m+1) - bmoinsun)
                cm.append(x)
        else:
            for j in range(1,k+1):
                cm.append((sum(PascalRow[i]*Rm(lesgammasprime[i]) *
                               Rm(touslescoeffs[m-i][j])
                               for i in range(1,m+1))
                           + cm[-1]
                           + sum(PascalRow[i]*Rm(lespuissancesdedprime[i]) *
                                 Rm(touslescoeffs[m-i][j-1])
                                 for i in range(1,m+1))
                           )
                          / Rm(b**(m+1) - bmoinsun)
                          )

        touslescoeffs.append(cm.copy())

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))

    # calcul parallèle des beta (sommes d'inverses de puissances)
    if showtimes:
        print(f"Calcul parallélisé des beta(m+1) avec maxworkers={maxworkers}")
        # We want to display some visual sign or progress
        # Find the largest multiple of maxworkers at most 1000
        # HOWEVER THIS COSTS IN CODE SIMPLICITY AND EFFICIENCY
        # ATTENTION maxworkers must be < 1000...
        mSize = (1000 // maxworkers) * maxworkers
        q, r = divmod(mSize, maxworkers)
        I = 0
        indices = []
        for i in range(r):
            indices.append(I)
            I += q + 1
        for i in range(maxworkers - r):
            indices.append(I)
            I += q
        # assert I == mSize, "check your math"
        indices.append(I)
        def map_para_beta(j):
            print(f"({j} occ.)...", end = ' ', flush = True)
            starttime = time.time()
            mend = 1
            L = [0]
            for i in range(Mmax//mSize):
                mbegin = 1 + i*mSize
                mend   = mbegin + mSize
                # TODO: find a way to not "instantiate" m range and keep
                #       it as a generator but I have to become more
                #       knowledgeable in Python and impact will be minor
                #       anyhow
                mrange = list(range(mbegin, mend))
                inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblockshifted[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(para_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{mend}", end = " ", flush= True)

            if mend < Mmax+1:
                mrange = list(range(mend, Mmax + 1))
                qlast, rlast = divmod(len(mrange), maxworkers)
                I = 0
                indiceslast = []
                for i in range(rlast):
                    indiceslast.append(I)
                    I += qlast + 1
                for i in range(maxworkers - rlast):
                    indiceslast.append(I)
                    I += qlast
                indiceslast.append(I)
                inputblocks = [(mrange[indiceslast[i]:indiceslast[i+1]], IndexToR,
                                maxblockshifted[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(para_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{Mmax+1} (done)", end = " ", flush=True)
            stoptime = time.time()
            print("{:.3f}s".format(stoptime - starttime))
            return L
    else:
        q, r = divmod(Mmax, maxworkers)
        I = 0
        indices = []
        for i in range(r):
            indices.append(I)
            I += q + 1
        for i in range(maxworkers - r):
            indices.append(I)
            I += q
        indices.append(I)
        def map_para_beta(j):
            L = [0]
            mrange = list(range(1, Mmax + 1))
            inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblockshifted[j])
                           for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(para_beta(inputblocks)))]
            L.extend(sum([result for result in results_1], []))
            return L

    lesbetas_maxblockshifted0 = map_para_beta(0)

    if k >= 1:
        lesbetas_maxblockshifted1 = map_para_beta(1)

    if k >= 2:
        lesbetas_maxblockshifted2 = map_para_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblockshifted3 = map_para_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblockshifted4 = map_para_beta(4)

    # boucle pour évaluer si all = True également les j < k
    Sk = []

    for j in range(0 if all else k, k+1):
        if showtimes:
            print(f"Calcul de l'approximation principale avec k={j}...",
                  end = ' ', flush = True)
            starttime = time.time()

        # calcul de la série positive de Burnol

        S = 0

        if j == 0:
            S = sum(1/R(x) for x in block1[0])
        elif j == 1:
            if d != 0:
                S = 1/R(d)

        if verbose:
            print("Somme du niveau 1 pour d = %s et j = %s:" % (d, j))
            print(S)

        if 2 < level:
            if j <= 2:
                S += sum(1/R(x) for x in block2[j])
            if verbose:
                print("Somme avec niveau 2 pour d = %s et j = %s:" % (d, j))
                print(S)

        if 3 < level:
            if j <= 3:
                S += sum(1/R(x) for x in block3[j])
            if verbose:
                print("Somme avec niveau 3 pour d = %s et j = %s:" % (d, j))
                print(S)

        # Attention à emploi de maxblockshifted ici
        S += b * (sum(sum(1/R(x) for x in maxblockshifted[i])
                      for i in range(1 + min(j,level))))

        if showtimes:
            stoptime = time.time()
            print("{:.3f}s".format(stoptime-starttime))

        if verbose:
            print(f"Somme ajustée de niveau {level} "
                  "avant incorporation de la série:")
            print(S)

        if verbose:
            print(f"On va utiliser {Mmax} termes de la série positive")

        if showtimes:
            print(f"Calcul de la série pour k={j}...", end = ' ', flush = True)
            starttime = time.time()

        # WE START WITH THE SMALLEST TERM CONTRIBUTING TO THE SERIES
        # Its sign will be set later.

        # Attention à emploi de maxblockshifted ici
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
                u, E = _shorten_small_real(lastterm)
                print(f"The {Mmax}th term is about {u:f} 10^{E} i.e. ",
                      end = "", flush = True)
            else:
                print(f"The {Mmax}th term is about {lastterm:.3e} i.e. ",
                      end = "", flush = True)

        # COMPUTATION OF THE MAIN SERIES BUILDING UP FROM SMALLEST TERMS
        # Rm will be the RealField. When m decreases Rm changes from time to
        # time regularly and automatically to use more bits.

        # Attention à emploi de maxblockshifted ici
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
            stoptime = time.time()
            print("{:.3f}s".format(stoptime-starttime))

        # NOW COMPUTE FINAL RESULT
        # This will later be trimmed from extra digits kept.
        S = S + bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = _shorten_small_real(ratio)
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

General information is also available in the irwin_v4_docstring
variable.

The variable "maxworkers" sets ncpus for @parallel usage. You
can set it and then reload.
{maxworkersinfostring}
"""
          )
