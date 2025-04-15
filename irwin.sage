# -*- mode: python ; coding: utf-8; -*-

"""Computes Irwin sums in high precision

Documentation is mostly in French.

Copyright (C) 2024 Jean-François Burnol
License: CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/

This file implements the formulas from the article 

    Measures for the summation of Irwin series
    J.-F. Burnol, February 2024

Change history:

v1: initial version. February 11, 2024

v2: February 26, 2024

    - make nbdigits an optional parameter with default value 34
    - modify a bit what gets printed under verbose=True option:
      only print a few significant digits of the "last term"
    - partially rewrite docstrings
    - add ``nbguardbits`` global variable to do computation with
      additional bits, defaults to 8 which is overly cautious.

v3: February 27, 2024

    - fix unfortunate stupid addition for incomprehensible reasons
      at v2 of

          assert type(b)==type(1)

      which broke usage of ``irwin(b,d,k)`` with a ``b`` if for example
      generated via range().

From a sage interactive session use

    load("irwin.sage")

"""

__version__ = "1.0.0"
__author__ = "Jean-François Burnol"
__copyright__ = "Copyright (C) 2024 Jean-François Burnol"
__license__ = "CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/"

nbguardbits = 8

def shorten_small_real(rr):
    """Clumsy way probably to display only magnitude order of a tiny real number

    """
    # I have to work around better this problem with float conversion
    # possibly giving zero due to low exponent
    x = RealField(53)(rr)
    s, _, _ = x.sign_mantissa_exponent()
    y = abs(x).log10()
    N = y.trunc()
    return s*10**float(y-N+1), N-1

print("Pré-calcul des coefficients du binôme jusque n=1000")

PascalC = [[1],[1,1]]

PascalC_mmax = 1000

for j in range(2,PascalC_mmax + 1):
    L = [1]
    for i in range(1,j):
        L.append(PascalC[j-1][i]+PascalC[j-1][i-1])
    L.append(1)
    PascalC.append(L)


def irwin(b, d, k, nbdigits=34, all=False, level=2, verbose=False, Mmax = -1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise algorithme de Burnol, série alternée de niveau 2, 3 ou 4

    :param int b: the integer base
    :param int d: the digit
    :param int k: the number of occurrences
    :param int nbdigits: (optional, default 34) target precision
    :param bool all: (optional, default ``False``
        if ``True``  print all irwin sums for ``j`` occurrences
        with ``j`` from ``0`` to ``k``.
        également pour tous les entiers inférieurs ou égaux à k.
    :param int level: (optional, default 2) level must be 2, 3 or 4
        le niveau souhaitée, 2, 3 ou 4 
    :param bool verbose: (optional, default ``False``)
        whether to print extra info
    :param int Mmax: (optional, default ``-1``)
        Forces to use that many terms of the Burnol series
    :rtype: :class:`sage.rings.real_mpfr.RealNumber`
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b

    Le meilleur choix entre level=2 (défaut) et level=3 dépend de b,
    nbdigits et k:

    - pour b=10, level=2 semble plus favorable jusqu'à
      environ nbdigits=1200/k.  Au-delà, c'est level=3 qui est conseillé.

    - pour de plus petites bases, level=3 voire même level=4 sont
      préférables à level=2 même pour nbdigits assez petit.

    Example:
    --------

    sage: irwin(10,9,3,all=True)
    (k= 0) 22.92067661926415034816365709437593
    (k= 1) 23.04428708074784831967594930973617
    (k= 2) 23.02604026596124378845022249787272
    (k= 3) 23.02585299837244431714290384468012
    23.02585299837244431714290384468012
    """

    global PascalC_mmax, PascalC  # I know this is bad

    assert 1 < level <= 4, "Le niveau (level) doit être 2 ou 3 ou 4"

    # SURTOUT NE PAS FAIRE if type(b) == type(1) !!!!
    # ça marche en pour des inputs directs mais pas pour "for b in range(B)"
    # J'avais bêtement ajouté cette erreur le samedi 24 février v2 sur arXiv
    assert b > 1, "%s doit être au moins 2" % b

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    nbbits = ceil((nbdigits+1)*log(10,2)) + nbguardbits
    R = RealField(nbbits)

    if verbose:
        print("Computations will be done with %s bits of precision" % nbbits)
        print("This includes nbguardbits (=%s) guard bits, which can be set" % nbguardbits)

    # NOTE: POUR k=0 ET d=1 CECI EST UNE GROSSE SURESTIMATION DU NOMBRE
    # DE TERMES RÉELLEMENT NÉCESSAIRES

    if Mmax == -1:
        Mmax = ceil(nbdigits/(level-1)/log(b,10))

    A = [i for i in range(b)]
    A.remove(d)

    N = b - 1

    lowblocks = [[0]]

    # attention dans le code originel pour les floats de Python
    # ce A1 avait déjà une conversion en float()
    A1 = [a for a in A if a !=0]

    # ATTENTION ! le lowblock pour 1 chiffre est défini comme pour k=0
    # il faudra rajouter le d pour les 1 <= j <= k plus bas
    block1 = []
    block1.append(A1)  # block1[0]: pas le chiffre d

    if d == 0:
        block1.append([])
    else:
        block1.append([d])

    block2 = []
    # block2[0] = nombres à deux chiffres sans le chiffre d
    # block2[1] = nombres à deux chiffres avec une occurrence de d
    # block2[2] = nombres à deux chiffres avec deux occurrences de d
    block2.append([b * x + a for x in block1[0] for a in A])  # k=0
    L = [b * x + d for x in block1[0]]
    L.extend([b * x + a for x in block1[1] for a in A])
    block2.append(L)  # k = 1
    if d == 0:
        block2.append([])  # k = 2
    else:
        block2.append([(b + 1) * d])  # k = 2

    # if verbose:
    #     print("cardinalité au niveau 2 pour k=0, 1, 2 (avec chiffre %s): %s, %s, %s" %
    #           (d, len(block2[0]), len(block2[1]), len(block2[2])))
        
    assert len(block2[0]) + len(block2[1]) + len(block2[2]) == b*(b-1), "Fix block2 bug!"

    if level == 2:
        maxblock = block2
    else:
        block3 = []
        block3.append([b * x + a for x in block2[0] for a in A])  # k=0
        L = [b * x + d for x in block2[0]]
        L.extend([b * x + a for x in block2[1] for a in A])
        block3.append(L)  # k=1
        L = [b * x + d for x in block2[1]]
        L.extend([b * x + a for x in block2[2] for a in A])
        block3.append(L)  # k=2
        if d == 0:
            block3.append([])
        else:
            block3.append([(b*b + b + 1) * d])  # k = 3

        # if verbose:
        #     print("cardinalité au niveau 3 pour k=0, 1, 2, 3 "
        #           "(avec chiffre %s): %s, %s, %s, %s" %
        #           (d,
        #            len(block3[0]),
        #            len(block3[1]),
        #            len(block3[2]),
        #            len(block3[3])
        #            )
        #           )

        assert (len(block3[0]) 
                + len(block3[1])
                + len(block3[2])
                + len(block3[3])) == b * b * (b-1), "Fix block3 bug!"

        if level == 3:
            maxblock = block3
        elif level == 4:
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

            maxblock = block4

            # if verbose:
            #     print("cardinalité au niveau 4 pour k=0, 1, 2, 3, 4 "
            #           "(avec chiffre %s): %s, %s, %s, %s, %s" %
            #           (d,
            #            len(block4[0]),
            #            len(block4[1]),
            #            len(block4[2]),
            #            len(block4[3]),
            #            len(block4[4])
            #        )
            #       )

            assert (len(block4[0]) 
                    + len(block4[1])
                    + len(block4[2])
                    + len(block4[3])
                    + len(block4[4])) == b * b * b * (b-1), "Fix block4 bug!"


    # Hésitation si en grande précision mieux de faire R(a**j) ou R(a)**j
    # ou de laisser la somme se faire avec des entiers

    lesgammas = [ N ]
    for j in range(1, Mmax+1):
        lesgammas.append(R(sum(a**j for a in A1)))
    
    lespuissancesded = [ 1 ]
    for j in range(1, Mmax+1):
        lespuissancesded.append(R(d**j))
    
    for j in range(PascalC_mmax+1,Mmax+1):
        L = [1]
        for i in range(1,j):
            L.append(PascalC[j-1][i]+PascalC[j-1][i-1])
        L.append(1)
        PascalC.append(L)

    PascalC_mmax = max(PascalC_mmax,Mmax)

    # calcul récursif des moments
    touslescoeffs = []

    # k = 0
    # hésitation sur est-ce que je dois convertir b en un R(b) pour les b**(m+1)
    # Rb = R(b)
    coeffs = [ R(b) ]
    for m in range(1, Mmax+1):
        coeffs.append(sum(PascalC[m][j]*lesgammas[j]*coeffs[m-j] for j in range(1,m+1))
                      / R(b**(m+1) - N)  # N = b-1
                      )  

    touslescoeffs.append(coeffs)

    for j in range(1,k+1):
        prevcoeffs = coeffs
        coeffs = [ R(b) ]
        for m in range(1, Mmax+1):
            coeffs.append((sum(PascalC[m][j]*lesgammas[j]*coeffs[m-j] for j in range(1,m+1))
                           + sum(PascalC[m][j]*lespuissancesded[j]*prevcoeffs[m-j] for j in range(0,m+1))
                           )
                          / R(b**(m+1)-N)  # N = b-1
                          )
        touslescoeffs.append(coeffs)
   
    # boucle pour évaluer si all = True également les j < k

    for j in range(0 if all else k, k+1):

        # calcul de la série alternée de Burnol

        S = 0

        if j == 0:
            S = sum(1/R(x) for x in A1)
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
        
        # jusqu'à j répétitions ; si j >= level, on s'arrête à level répétitions max
        S += b * (sum(sum(1/R(x) for x in maxblock[i]) for i in range(1 + min(j,level))))

        if verbose:
            print("Somme ajustée de niveau %s avant corrections géométriques:" %
                  level)
            print(S)

        if verbose:
            print("On va utiliser %d termes de la série alternée" % Mmax)

        bubu = touslescoeffs[j][Mmax] * sum(1/R(n ** (Mmax+1)) for n in maxblock[0])
        if j >= 1: 
            bubu += touslescoeffs[j-1][Mmax] * sum(1/R(n ** (Mmax+1)) for n in maxblock[1])
        if j >= 2:
            bubu += touslescoeffs[j-2][Mmax] * sum(1/R(n ** (Mmax+1)) for n in maxblock[2])
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[j-3][Mmax] * sum(1/R(n ** (Mmax+1)) for n in maxblock[3])
        if (level == 4) and (j >= 4):
            bubu += touslescoeffs[j-4][Mmax] * sum(1/R(n ** (Mmax+1)) for n in maxblock[4])
                
        if verbose:
            lastterm = -bubu if Mmax&1 else bubu
            if float(lastterm) == 0.:
                u, E = shorten_small_real(lastterm)
                print("The %sth term is about %f times 10^%s and it represents" % (Mmax, u, E))
            else:
                print("The %sth term is about %.3e and it represents" % (Mmax, lastterm))

        for m in range(Mmax-1, 0, -1):  # last one is m=1
            bubu = -bubu
            bubu += touslescoeffs[j][m] * sum(1/R(n ** (m+1)) for n in maxblock[0])
            if j >= 1:
                bubu += touslescoeffs[j-1][m] * sum(1/R(n ** (m+1)) for n in maxblock[1])
            if j >= 2:
                bubu += touslescoeffs[j-2][m] * sum(1/R(n ** (m+1)) for n in maxblock[2])
            if (level > 2) and (j >= 3):
                bubu += touslescoeffs[j-3][m] * sum(1/R(n ** (m+1)) for n in maxblock[3])
            if (level == 4) and (j >= 4):
                bubu += touslescoeffs[j-4][m] * sum(1/R(n ** (m+1)) for n in maxblock[4])

        S = S - bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = shorten_small_real(ratio)
                print("a fraction %.3f times 10^%s of complete sum." % (u, E))
            else:
                print("a fraction %.3e of complete sum." % ratio)

            print("La somme de m=1 à %s vaut" % Mmax)
            print(-bubu)

        if all:
            print("(k=%2d) %s" % (j, S))

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    Rnbdigits = RealField(nbbits-nbguardbits)
    return Rnbdigits(S)


def irwinpos(b, d, k, nbdigits=34, all=False, level=2, verbose=False, Mmax = -1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise algorithme de Burnol, série positive de niveau 2, 3 ou 4

    :param int b: the integer base
    :param int d: the digit
    :param int k: the number of occurrences
    :param int nbdigits: (optional, default 34) target precision
    :param bool all: (optional, default ``False``
        if ``True``  print all irwin sums for ``j`` occurrences
        with ``j`` from ``0`` to ``k``.
        également pour tous les entiers inférieurs ou égaux à k.
    :param int level: (optional, default 2) level must be 2, 3 or 4
        le niveau souhaitée, 2, 3 ou 4 
    :param bool verbose: (optional, default ``False``)
        whether to print extra info
    :param int Mmax: (optional, default ``-1``)
        Forces to use that many terms of the Burnol series
    :rtype: :class:`sage.rings.real_mpfr.RealNumber`
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b

    Le meilleur choix entre level=2 (défaut) et level=3 dépend de b,
    nbdigits et k:

    - pour b=10, level=2 semble plus favorable jusqu'à
      environ nbdigits=1200/k.  Au-delà, c'est level=3 qui est conseillé.

    - pour de plus petites bases, level=3 voire même level=4 sont
      préférables à level=2 même pour nbdigits assez petit.

    Example:
    --------

    sage: irwinpos(10,9,3,all=True)
    (k= 0) 22.92067661926415034816365709437593
    (k= 1) 23.04428708074784831967594930973618
    (k= 2) 23.02604026596124378845022249787272
    (k= 3) 23.02585299837244431714290384468012
    23.02585299837244431714290384468012
    """

    global PascalC_mmax, PascalC  # I know this is bad

    assert 1 < level <= 4, "Le niveau (level) doit être 2, 3 ou 4"

    # SURTOUT NE PAS FAIRE if type(b) == type(1) !!!!
    # ça marche en pour des inputs directs mais pas pour "for b in range(B)"
    # J'avais bêtement ajouté cette erreur le samedi 24 février v2 sur arXiv
    assert b > 1, "%s doit être au moins 2" % b

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    nbbits = ceil((nbdigits+1)*log(10,2)) + nbguardbits
    R = RealField(nbbits)

    if verbose:
        print("Computations will be done with %s bits of precision" % nbbits)
        print("This includes nbguardbits (=%s) guard bits, which can be set" % nbguardbits)

    # NOTE: POUR k=0 ET d=1 CECI EST UNE GROSSE SURESTIMATION DU NOMBRE
    # DE TERMES RÉELLEMENT NÉCESSAIRES

    if Mmax == -1:
        Mmax = ceil(nbdigits/(level-1)/log(b,10))

    A = [i for i in range(b)]
    A.remove(d)

    N = b- 1

    lowblocks = [[0]]

    # attention dans le code originel pour les floats de Python
    # ce A1 avait déjà une conversion en float()
    A1 = [a for a in A if a !=0]

    # ATTENTION ! le lowblock pour 1 chiffre est défini comme pour k=0
    # il faudra rajouter le d pour les 1 <= j <= k plus bas
    block1 = []
    block1.append(A1)  # block1[0]: pas le chiffre d

    if d == 0:
        block1.append([])
    else:
        block1.append([d])

    block2 = []
    # block2[0] = nombres à deux chiffres sans le chiffre d
    # block2[1] = nombres à deux chiffres avec une occurrence de d
    # block2[2] = nombres à deux chiffres avec deux occurrences de d
    block2.append([b * x + a for x in block1[0] for a in A])  # k=0
    L = [b * x + d for x in block1[0]]
    L.extend([b * x + a for x in block1[1] for a in A])
    block2.append(L)  # k = 1
    if d == 0:
        block2.append([])  # k = 2
    else:
        block2.append([(b + 1) * d])  # k = 2

    # if verbose:
    #     print("cardinalité au niveau 2 pour k=0, 1, 2 (avec chiffre %s): %s, %s, %s" %
    #           (d, len(block2[0]), len(block2[1]), len(block2[2])))
        
    assert len(block2[0]) + len(block2[1]) + len(block2[2]) == b*(b-1), "Fix block2 bug!"

    if level == 2:
        maxblock = block2
    else:
        block3 = []
        block3.append([b * x + a for x in block2[0] for a in A])  # k=0
        L = [b * x + d for x in block2[0]]
        L.extend([b * x + a for x in block2[1] for a in A])
        block3.append(L)  # k=1
        L = [b * x + d for x in block2[1]]
        L.extend([b * x + a for x in block2[2] for a in A])
        block3.append(L)  # k=2
        if d == 0:
            block3.append([])
        else:
            block3.append([(b*b + b + 1) * d])  # k = 3

        # if verbose:
        #     print("cardinalité au niveau 3 pour k=0, 1, 2, 3 "
        #           "(avec chiffre %s): %s, %s, %s, %s" %
        #           (d,
        #            len(block3[0]),
        #            len(block3[1]),
        #            len(block3[2]),
        #            len(block3[3])
        #            )
        #           )

        assert (len(block3[0]) 
                + len(block3[1])
                + len(block3[2])
                + len(block3[3])) == b * b * (b-1), "Fix block3 bug!"

        if level == 3:
            maxblock = block3
        elif level == 4:
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

            maxblock = block4

            # if verbose:
            #     print("cardinalité au niveau 4 pour k=0, 1, 2, 3, 4 "
            #           "(avec chiffre %s): %s, %s, %s, %s, %s" %
            #           (d,
            #            len(block4[0]),
            #            len(block4[1]),
            #            len(block4[2]),
            #            len(block4[3]),
            #            len(block4[4])
            #        )
            #       )

            assert (len(block4[0]) 
                    + len(block4[1])
                    + len(block4[2])
                    + len(block4[3])
                    + len(block4[4])) == b * b * b * (b-1), "Fix block4 bug!"


    # ATTENTION que la série positive a des récurrences avec b-1-d à la place de d
    Aprime = [i for i in range(b)]
    dprime = b-1-d
    Aprime.remove(dprime)

    lesgammasprime = [N]
    for j in range(1, Mmax+1):
        lesgammasprime.append(R(sum(a**j for a in Aprime if a != 0)))

    Rdprime=R(dprime)    
    lespuissancesdedprime = [ 1 ]
    for j in range(1, Mmax+1):
        lespuissancesdedprime.append(R(dprime**j))
    
    for j in range(PascalC_mmax+1,Mmax+1):
        L = [1]
        for i in range(1,j):
            L.append(PascalC[j-1][i]+PascalC[j-1][i-1])
        L.append(1)
        PascalC.append(L)

    PascalC_mmax = max(PascalC_mmax, Mmax)

    # calcul récursif des moments
    touslescoeffs = []

    # k = 0
    # hésitation sur utilisation de Rb pour les b**(m+1)
    # Rb = R(b)
    coeffs = [ R(b) ]
    for m in range(1, Mmax+1):
        coeffs.append((R(b**(m+1))  # terme additionnel pour la récurrence des v
                       + sum(PascalC[m][j]*lesgammasprime[j]*coeffs[m-j] for j in range(1,m+1))
                       )/ R(b**(m+1) - N)  # N = b-1
                    )

    touslescoeffs.append(coeffs)

    for j in range(1,k+1):
        prevcoeffs = coeffs
        coeffs = [ R(b) ]
        for m in range(1, Mmax+1):
            coeffs.append((sum(PascalC[m][j]*lesgammasprime[j]*coeffs[m-j] for j in range(1,m+1))
                           + sum(PascalC[m][j]*lespuissancesdedprime[j]*prevcoeffs[m-j] for j in range(0,m+1))
                           )
                          / R(b**(m+1) - N)  # N = b-1
                          )
        touslescoeffs.append(coeffs)
   
    # boucle pour évaluer si all = True également les j < k

    for j in range(0 if all else k, k+1):

        # calcul de la série positive de Burnol

        S = 0

        if j == 0:
            S = sum(1/R(x) for x in A1)
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
               
        S += b * (sum(sum(1/R(x+1) for x in maxblock[i]) for i in range(1 + min(j,level))))

        if verbose:
            print("Somme ajustée de niveau %s avant corrections géométriques:" %
                  level)
            print(S)

        if verbose:
            print("On va utiliser %d termes de la série alternée" % Mmax)

        bubu = touslescoeffs[j][Mmax] * sum(1/R((n + 1) ** (Mmax+1)) for n in maxblock[0])
        if j >= 1: 
            bubu += touslescoeffs[j-1][Mmax] * sum(1/R((n + 1) ** (Mmax+1)) for n in maxblock[1])
        if j >= 2:
            bubu += touslescoeffs[j-2][Mmax] * sum(1/R((n + 1) ** (Mmax+1)) for n in maxblock[2])
        if (level > 2) and (j >= 3):
            bubu += touslescoeffs[j-3][Mmax] * sum(1/R((n + 1) ** (Mmax+1)) for n in maxblock[3])
        if (level == 4) and (j >= 4):
            bubu += touslescoeffs[j-4][Mmax] * sum(1/R((n + 1) ** (Mmax+1)) for n in maxblock[4])
                
        if verbose:
            lastterm = -bubu if Mmax&1 else bubu
            if float(lastterm) == 0.:
                u, E = shorten_small_real(lastterm)
                print("The %sth term is about %f times 10^%s and it represents" % (Mmax, u, E))
            else:
                print("The %sth term is about %.3e and it represents" % (Mmax, lastterm))

        for m in range(Mmax-1, 0, -1):  # last one is m=1
            bubu += touslescoeffs[j][m] * sum(1/R((n + 1) ** (m+1)) for n in maxblock[0])
            if j >= 1:
                bubu += touslescoeffs[j-1][m] * sum(1/(R(n + 1) ** (m+1)) for n in maxblock[1])
            if j >= 2:
                bubu += touslescoeffs[j-2][m] * sum(1/(R(n + 1) ** (m+1)) for n in maxblock[2])
            if (level > 2) and (j >= 3):
                bubu += touslescoeffs[j-3][m] * sum(1/R((n + 1) ** (m+1)) for n in maxblock[3])
            if (level == 4) and (j >= 4):
                bubu += touslescoeffs[j-4][m] * sum(1/R((n + 1) ** (m+1)) for n in maxblock[4])

        S = S + bubu

        if verbose:
            ratio = lastterm/S
            if float(ratio) == 0.:
                u, E = shorten_small_real(ratio)
                print("a fraction %.3f times 10^%s of complete sum." % (u, E))
            else:
                print("a fraction %.3e of complete sum." % ratio)

            print("La somme de m=1 à %s vaut" % Mmax)
            print(-bubu)

        if all:
            print("(k=%2d) %s" % (j, S))

    if verbose:
        print("b = %s, d = %s, k = %s, level = %s" % (b, d, k, level))

    Rnbdigits = RealField(nbbits-nbguardbits)
    return Rnbdigits(S)


def irwintest():
    """Comparaison avec données qu'on peut récupérer dans article de Baillie, 2008

    La plupart en base 10 et avec 20 chiffre après virgule.

    Inclut aussi un test de comparaison avec des résultats à 100 chiffres
    calculés en 2012 pour Kempner.  Comme ceux-ci étaient tronqués, pas arrondis,
    on a arrondi à 99 chiffres après la virgule.

    Teste les niveaux 2, 3 et 4.  Ce dernier ralentit significativement le test.
    """
    bailliedata = [
        [ 23.10344790942054161603, 23.02673534156912696109, 23.02586068273551997642 ],
        [ 16.17696952812344426658, 23.16401859427283204085, 23.02727628635600571224 ],
        [ 19.25735653280807222453, 23.08826066275634239334, 23.02648597376847065598 ],
        [ 20.56987795096123037108, 23.06741088193023010242, 23.02627319066793505960 ],
        [ 21.32746579959003668664, 23.05799241338182439576, 23.02617788539260017317 ],
        [ 21.83460081229691816341, 23.05272889453011749904, 23.02612487531564760861 ],
        [ 22.20559815955609188417, 23.04940997329550055704, 23.02609154986488712587 ],
        [ 22.49347531170594539818, 23.04714619019864185083, 23.02606886491441507436 ],
        [ 22.72636540267937060283, 23.04551390798215553342, 23.02605253084569367648 ],
        [ 22.92067661926415034816, 23.04428708074784831968, 23.02604026596124378845 ],
    ]

    R22 = RealField(22*log(10,2))
    for d in range(10):
        for k in range(3):
            y = R22(bailliedata[d][k])

            x = irwin(10,d,k,22); print(x, " (10,%s,%s,22)" % (d,k));
            assert abs(x-y) < 1e-20,\
                "irwin(10,%s,%s,22), delta=%s" % (d, k, float(x-y))

            x = irwin(10,d,k,22,level=3); print(x, " (level=3)");
            assert abs(x-y) < 1e-20,\
                "irwin(10,%s,%s,22,level=3), delta=%s" % (d, k, float(x-y))

            x = irwin(10,d,k,22,level=4); print(x, " (level=4)");
            assert abs(x-y) < 1e-19,\
                "irwin(10,%s,%s,22,level=4), delta=%s" % (d, k, float(x-y))

            x = irwinpos(10,d,k,22); print(x, " (pos)");
            assert abs(x-y) < 1e-20,\
                "irwinpos(10,%s,%s,22), delta=%s" % (d, k, float(x-y))

            x = irwinpos(10,d,k,22,level=3); print(x, " (pos, level=3)");
            assert abs(x-y) < 1e-20,\
                "irwinpos(10,%s,%s,22,level=3), delta=%s" % (d, k, float(x-y))

            x = irwinpos(10,d,k,22,level=4); print(x, " (pos, level=4)");
            assert abs(x-y) < 1e-19,\
                "irwinpos(10,%s,%s,22,level=4), delta=%s" % (d, k, float(x-y))

    def test(args, expected, eps):
        R = RealField((args[-1]+1)*log(10,2))

        x = irwin(*args); print(x, " ", args);
        delta = x - R(expected)
        assert abs(delta) < eps, "irwin(%s), delta=%s" % (args,delta) 

        x = irwinpos(*args); print(x, " (pos) ");
        delta = x - R(expected)
        #  print(float(delta))
        assert abs(delta) < eps, "irwinpos(%s), delta=%s" % (args,delta) 

        x = irwin(*args, level=3); print(x, " (level 3)");
        delta = x - R(expected)
        # surprenant pour le cas b=10,d=0,k=0,nbdigits=101 j'ai eu à un moment
        # un delta de -4.57e-100 mais pourtant l'affichage se termine avec
        # 6139338 pas 6139339 avec pour cette affichage un ulp de 1e-99
        # print(float(delta))
        assert abs(delta) < eps, "irwin(%s,level=3), delta=%s" % (args,delta) 

        x = irwinpos(*args, level=3); print(x, " (pos level 3)");
        delta = x - R(expected)
        assert abs(delta) < eps, "irwin(%s,level=3), delta=%s" % (args,delta) 

        # Au niveau 4 j'autorise un plus grand écart car les sommes avec
        # pour de l'ordre de b^4 termes accumulent des erreurs d'arrondis

        x = irwin(*args, level=4); print(x, " (level 4)");
        delta = x - R(expected)
        # print(float(delta))
        assert abs(delta) < 2*eps, "irwin(%s,level=4), delta=%s" % (args,delta) 

        x = irwinpos(*args, level=4); print(x, " (pos level 4)");
        delta = x - R(expected)
        # print(float(delta))
        assert abs(delta) < 2*eps, "irwin(%s,level=4), delta=%s" % (args,delta) 

    test((10,0,3,17),  23.025851037148538, 1e-15)
    test((10,0,10,17), 23.025850929940457, 1e-15)
    test((2,0,0,21),    1.60669515241529176378, 1e-20)
    test((2,0,1,26),    1.4625907350443646995461454, 1e-25)
    test((2,1,0,20), 0., 1e-30)
    test((2,1,1,20), 2., 1e-19)

    # les données de 2012 sont tronquées pas arrondies, donc ici j'ai arrondi
    # à 99 chiffres après la virgule
    for y in [
            (22.920676619264150348163657094375931914944762436998481568541998356572156338189911129445626037448201899, 9),
            (22.726365402679370602833644156742557889210702616360219843536376162400468201751348127010562165158922478, 8),
            (22.493475311705945398176226915339775974005915541672512361791460444071051200950740851432220823450021919, 7),
            (22.205598159556091884167380480007527105193856106668463270276938233053228350891247526347776997405891493, 6),
            (21.834600812296918163407235040609182717846567515013918291679359184250862668822938357772138319329254881, 5),
            (21.327465799590036686639401486939512843750951703270021817251189541977884272451335375381201302840693548, 4),
            (20.569877950961230371075217419053111414153869674730783489508528500267294996193803500590474940806035350, 3),
            (19.257356532808072224532776770194454115526053831154870149868362949104309016019551809280546221128442864, 2),
            (16.176969528123444266579603880364009305567219790763133864516906490836362988999996456388862146266850286, 1),
            (23.103447909420541616034054043325598138302800005282141886723094772738750796061419426359201910526139339, 0),
    ]:
        test((10, y[1], 0, 101), y[0], 1e-99)


if __name__ == "__main__":
    print("""
irwin(b,d,k,nbdigits): calcule avec N chiffres décimaux de précision
                la somme des 1/n pour n ayant le chiffre d
                en base b exactement k fois.  Le calcul utilise
                la série alternée de Burnol.
irwinpos(b,d,k,nbdigits): idem avec la série positive de Burnol.

Example of use:

    sage: irwin(10,9,10,52)
    23.02585092994045684021991777757086344589171431282125

The nbdigits parameter is optional and defaults to 34 decimal digits.
Issue help(irwin) or help(irwinpos) for additional info.
"""
          )


