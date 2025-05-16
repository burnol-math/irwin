# -*- mode: python ; coding: utf-8; -*-

# irwin_v4.sage OBSOLETE and stripped from most comments
# Use via load("irwin_v4.sage") in sage interactive mode

__version__  = "1.5.0"
__date__     = "2025/05/16"
__filename__ = "irwin_v4.sage"

irwin_v4_docstring = """
This file is obsolete.  Use irwin_v5.sage or later rather.

Thanks to Arnaud Bodin and Yusuf Emin Akpinar for interactions.

Copyright (C) 2025 Jean-François Burnol
License: CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/
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
def _v4_um(a, P, G, T, R, j, m):
    return sum(P[i]*R(G[i])*R(T[m-i][j]) for i in range(a,m+1,maxworkers))

def _v4_beta_aux(m, R, nblock):
    return sum(1/R(n ** (m+1)) for n in nblock)
@parallel(ncpus=maxworkers)
def _v4_beta(mlist, IR, nblock):
    return list(_v4_beta_aux(m, IR[m], nblock) for m in mlist)

nbguardbits = 12

def _v4_shorten_small_real(rr):
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


def _v4_setup_realfields(nbdigits, PrecStep, b, level, Mmax=-1):
    """Preparation of an array mapping index to a RealField
    """

    # Try to guarantee we will have nbdigits radix 10 digits correct.
    nbbits_final = ceil((nbdigits+1)*log(10,2))
    Rfinal = RealField(nbbits_final)
    if True:
        nbbits = nbbits_final + nbguardbits
    else:
        nbbits = 32 * ceil((nbbits_final + nbguardbits)/32)
    R = RealField(nbbits)

    _Mmax = floor((nbbits - nbguardbits/2)/(level-1)/log(b,2))
    NbOfPrec = 1 + floor((nbbits - nbguardbits/2)/PrecStep)
    LesReels = [R]
    for j in range(1, NbOfPrec):
        LesReels.append(RealField(nbbits - j * PrecStep))
    IndexToPrecRatio = PrecStep / (level - 1) / log(b,2)  # not an integer !
    IndexToR = []
    oldindexbound = 0
    for j in range(NbOfPrec):
        newindexbound = ceil((j+1)*IndexToPrecRatio)
        IndexToR.extend([LesReels[j]] * (newindexbound - oldindexbound))
        oldindexbound = newindexbound
    IndexToR[1] = R  # LesReels[0]
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


def _v4_setup_blocks(b, d, level):
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


def _v4_umtimeinfo(single, multi, para, wrkrs, m):
    """Auxiliary shared between irwin() and irwinpos()
    """
    if multi < single:
        if para:
            print(f"... poursuite car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle (m={m})")
        else:
            print(f"... basculement car {single:.3f}s>{multi:.3f}s"
                  f" en parallèle (maxworkers={wrkrs}, m={m})")
    else:
        if para:
            print(f"... on quitte car {single:.3f}s<{multi:.3f}s"
                  f" l'exécution parallèle (m={m}) ")
        else:
            print(f"... pas utile ({single:.3f}s<{multi:.3f}s)"
                  f" d'exécuter en parallèle (m={m}) ")

def irwin(b, d, k,
          nbdigits=34,
          level=3,
          PrecStep=500,
          all=False,
          showtimes=False,
          verbose=False,
          Mmax=-1):
    """Somme d'Irwin pour b, d, k avec nbdigits chiffres décimaux (en tout).

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
     NbOfPrec) = _v4_setup_realfields(nbdigits, PrecStep, b, level, Mmax)

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

    blocks = _v4_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level >2:
        block3 = blocks[2]
        if level >3:
            block4 = blocks[3]
    maxblock = blocks[-1]
    lesgammas = [ bmoinsun ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lesgammas.append(Rj(sum(a**j for a in block1[0])))

    lespuissancesded = [ 1 ]
    for j in range(1, Mmax+1):
        Rj = IndexToR[j]
        lespuissancesded.append(Rj(d**j))

    if showtimes:
        stoptime = time.time()
        print("{:.3f}s".format(stoptime - starttime))

    if showtimes:
        if k == 0:
            print(f"Calcul des u_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des u_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.time()

    touslescoeffs = [ [R(b)] * (k+1) ]
    c1 = [ lesgammas[1] * R(b) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(((lesgammas[1] + d) * R(b) + c1[-1])/R(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRow = [1, 1]
    useparallel = False
    if Mmax > 511 and showtimes and verbose:
        print("    L'algorithme teste tous les 512 coefficients s'il est bénéfique")
        print("    d'utiliser la procédure décorée par @parallel.  Ce test n'est fait")
        print("    que pour le calcul des u_{0;m}.  Le premier nombre indiqué sera")
        print("    la durée du calcul normal, le second avec @parallel, pour u_{0;m}.")

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
                results = list(_v4_um(((a, PascalRow, lesgammas, touslescoeffs, Rm, 0, m,)
                                       for a in range(1, maxworkers + 1))))
                partial_sums = [result[1] for result in results]
                cm = [ sum(partial_sums) / Rm(b**(m+1) - bmoinsun) ]
            else:
                cm = [ sum(PascalRow[i]*Rm(lesgammas[i])*Rm(touslescoeffs[m-i][0])
                           for i in range(1,m+1)) / Rm(b**(m+1) - bmoinsun) ]
        else:
            localstarttime = time.time()
            results = list(_v4_um(((a, PascalRow, lesgammas, touslescoeffs, Rm, 0, m,)
                                   for a in range(1, maxworkers + 1))))
            partial_sums = [result[1] for result in results]
            x = [ sum(partial_sums) / Rm(b**(m+1) - bmoinsun) ]
            multitime = time.time() - localstarttime

            localstarttime = time.time()
            cm = [ sum(PascalRow[i]*Rm(lesgammas[i])*Rm(touslescoeffs[m-i][0])
                       for i in range(1,m+1)) / Rm(b**(m+1) - bmoinsun) ]

            singletime = time.time() - localstarttime

            if showtimes:
                _v4_umtimeinfo(singletime, multitime,
                               useparallel, maxworkers, m)

            useparallel = (multitime < singletime)

        if useparallel:
            for j in range(1,k+1):
                results = list(_v4_um(((a, PascalRow, lesgammas,
                                        touslescoeffs, Rm, j, m,)
                                       for a in range(1, maxworkers + 1))))
                x = sum(p[1] for p in results)
                x += cm[-1]
                results = list(_v4_um(((a, PascalRow, lespuissancesded,
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
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} (fait) "
              + f"{stoptime-starttime:.3f}s")

    # calcul parallèle des beta (sommes d'inverses de puissances)
    if showtimes:
        print("Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
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
        def map__v4_beta(j):
            print(f"... ({j} occ.) ", end = "", flush = True)
            starttime = time.time()
            mbegin = 1
            mend = 1
            L = [0]
            for rep in range(Mmax//mSize):
                mend = mbegin + mSize
                mrange = list(range(mbegin, mend))
                inputblocks = [(mrange[indices[i]:indices[i+1]],
                                IndexToR,
                                maxblock[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(_v4_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{mend}", end = " ", flush= True)
                if (rep + 1) & 7 == 0:
                    print(f"\n" + " " * 12, end = " ")
                mbegin = mend

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
                             in sorted(list(_v4_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{Mmax+1} (fait)", end = " ", flush=True)
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
        def map__v4_beta(j):
            L = [0]
            mrange = list(range(1, Mmax + 1))
            inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblock[j])
                           for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v4_beta(inputblocks)))]
            L.extend(sum([result for result in results_1], []))
            return L

    lesbetas_maxblock0 = map__v4_beta(0)

    if k >= 1:
        lesbetas_maxblock1 = map__v4_beta(1)

    if k >= 2:
        lesbetas_maxblock2 = map__v4_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblock3 = map__v4_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblock4 = map__v4_beta(4)

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
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
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
                u, E = _v4_shorten_small_real(lastterm)
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
                u, E = _v4_shorten_small_real(ratio)
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


def irwinpos(b, d, k,
             nbdigits=34,
             level=3,
             PrecStep=500,
             all=False,
             showtimes=False,
             verbose=False,
             Mmax=-1):
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
        starttime = time.time()

    (nbbits,
     nbbits_final,
     R,
     Rfinal,
     Mmax,
     IndexToR,
     NbOfPrec) = _v4_setup_realfields(nbdigits, PrecStep, b, level, Mmax)

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

    blocks = _v4_setup_blocks(b, d, level)
    block1 = blocks[0]
    block2 = blocks[1]
    if level >2:
        block3 = blocks[2]
        if level >3:
            block4 = blocks[3]
    maxblock = blocks[-1]
    maxblockshifted = [[ n + 1  for n in L] for L in maxblock]

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
        if k == 0:
            print(f"Calcul des v_{{0;m}} pour m<={Mmax} ...")
        else:
            print(f"Calcul des v_{{j;m}} pour j<={k} et m<={Mmax} ...")
        starttime = time.time()

    touslescoeffs = [ [R(b)] * (k+1) ]
    # attention to b * b  extra in first one
    c1 = [ (b * b + lesgammasprime[1] * R(b)) / (b * b - bmoinsun) ]
    for j in range(1, k+1):
        c1.append(( (lesgammasprime[1] + dprime) * R(b) + c1[-1])/R(b * b - bmoinsun))
    touslescoeffs.append(c1)

    PascalRow = [1, 1]
    useparallel = False
    if Mmax > 511 and showtimes and verbose:
        print("    L'algorithme teste tous les 512 coefficients s'il est bénéfique")
        print("    d'utiliser la procédure décorée par @parallel.  Ce test n'est fait")
        print("    que pour le calcul des v_{0;m}.  Le premier nombre indiqué sera")
        print("    la durée du calcul normal, le second avec @parallel, pour v_{0;m}.")

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
                results = list(_v4_um(((a, PascalRow, lesgammasprime,
                                        touslescoeffs, Rm, 0, m,)
                                       for a in range(1, maxworkers + 1))))
                partial_sums = [result[1] for result in results]
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
            results = list(_v4_um(((a, PascalRow, lesgammasprime,
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
                _v4_umtimeinfo(singletime, multitime,
                               useparallel, maxworkers, m)

            useparallel = (multitime < singletime)

        if useparallel:
            for j in range(1,k+1):
                results = list(_v4_um(((a, PascalRow, lesgammasprime,
                                        touslescoeffs, Rm, j, m,)
                                       for a in range(1, maxworkers + 1))))
                x = sum(p[1] for p in results)
                x += cm[-1]
                results = list(_v4_um(((a, PascalRow, lespuissancesdedprime,
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
        print(f"... m<={Mmax}{f' et j<={k}' if k>0 else ''} (fait) "
              + f"{stoptime-starttime:.3f}s")

    # calcul parallèle des beta (sommes d'inverses de puissances)
    if showtimes:
        print(f"Calcul parallélisé des beta(m+1) avec "
              f"maxworkers={maxworkers} ...")
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
        def map__v4_beta(j):
            print(f"... ({j} occ.) ", end = "", flush = True)
            starttime = time.time()
            mbegin = 1
            mend = 1
            L = [0]
            for rep in range(Mmax//mSize):
                mend   = mbegin + mSize
                mrange = list(range(mbegin, mend))
                inputblocks = [(mrange[indices[i]:indices[i+1]],
                                IndexToR,
                                maxblockshifted[j])
                               for i in range(maxworkers)]
                results_1 = [result[1] for result
                             in sorted(list(_v4_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{mend}", end = " ", flush= True)
                if (rep + 1) & 7 == 0:
                    print(f"\n" + " " * 12, end = " ")
                mbegin = mend

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
                             in sorted(list(_v4_beta(inputblocks)))]
                L.extend(sum([result for result in results_1], []))
                print(f"m<{Mmax+1} (fait)", end = " ", flush=True)
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
        def map__v4_beta(j):
            L = [0]
            mrange = list(range(1, Mmax + 1))
            inputblocks = [(mrange[indices[i]:indices[i+1]], IndexToR, maxblockshifted[j])
                           for i in range(maxworkers)]
            results_1 = [result[1] for result
                         in sorted(list(_v4_beta(inputblocks)))]
            L.extend(sum([result for result in results_1], []))
            return L

    lesbetas_maxblockshifted0 = map__v4_beta(0)

    if k >= 1:
        lesbetas_maxblockshifted1 = map__v4_beta(1)

    if k >= 2:
        lesbetas_maxblockshifted2 = map__v4_beta(2)

    if (k >= 3) and (level > 2):
        lesbetas_maxblockshifted3 = map__v4_beta(3)

    if (k >= 4) and (level > 3):
        lesbetas_maxblockshifted4 = map__v4_beta(4)

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
            print("\nSomme du niveau 1 pour d = %s et j = %s:" % (d, j))
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
                u, E = _v4_shorten_small_real(lastterm)
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
                u, E = _v4_shorten_small_real(ratio)
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
    print(f"""\
WARNING: this 'v4' version is OBSOLETE.  Use 'v5' or later rather.
WARNING: this 'v4' version is OBSOLETE.  Use 'v5' or later rather.
WARNING: this 'v4' version is OBSOLETE.  Use 'v5' or later rather.
WARNING: well, you probably got it by now.
"""
          )
