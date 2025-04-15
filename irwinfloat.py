# -*- mode: python ; coding: utf-8; -*-

"""This module computes Irwin sums using Python floats

Documentation is in French.

date: vendredi 9 février 2024
authors: Jean-François Burnol

Copyright (C) 2024 Jean-François Burnol
License: CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/

Ce fichier définit deux procédures irwinfloat() et irwinposfloat() qui
implémentent avec des floats en Python les deux algorithmes de Burnol pour le
calcul des sommes de Kempner-Irwin, en toute base.

Exécuter en ligne de commande

    python -i irwinfloat.py

et suivre les instructions... (voir au bas du fichier)

Pour des calculs avec plus de chiffres utiliser Sagemath avec irwin.sage

"""

__version__ = "1.0.0"
__author__ = "Jean-François Burnol"
__copyright__ = "Copyright (C) 2024 Jean-François Burnol"
__license__ = "CC BY-SA 4.0 https://creativecommons.org/licenses/by-sa/4.0/"

from math import log10, ceil

# pour les coefficients du binôme, mais certains collègues n'ont pas scipy
# from scipy.special import binom
# https://stackoverflow.com/a/31435571

# calcul des coefficients du binôme
# on n'en a besoin que jusque Pascal(51,j) au maximum.
# 247959266474052 = binom(51,25) a 15 chiffres, calculer
# en floats n'apporterait aucun gain de vitesse
# Pour algorithme positif on prend un terme de plus par sécurité
# donc jusque Pascal(52,26)
Pascal = [[1],[1,1]]

for j in range(2,55):
    L = [1]
    for i in range(1,j):
        L.append(Pascal[j-1][i]+Pascal[j-1][i-1])
    L.append(1)
    Pascal.append(L)


def irwinfloat(b, d, k, all=False, verbose=False, Mmax = -1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise algorithme de Burnol, série alternée de niveau 2

    :param b: int
        la base
    :param d: int
        le chiffre
    :param k: int
        l'entier k
    :param all: bool,optional (False)
        si True, on calcule et on affiche
        également pour tous les entiers inférieurs ou égaux à k.
    :param verbose: bool,optional (False)
        si True affiche informations supplémentaires
    :param Mmax: int,optional (-1)
        si positif ou nul, donne le nombre de termes de la série de Burnol
        à utiliser
    :rtype: float
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b

    """

    assert b>1, "%s doit être au moins 2" % b

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    if Mmax == -1:
        Mmax = ceil(15/log10(b)) +1  # +1 de sécurité utile que pour des grands b
                                     # mais pour ceux-là de toute façon algo pas praticable

    if verbose:
        print("Mmax = ", Mmax)

    # # auxiliaire
    # def fb(a):
    #     return lambda x : b * x + a

    A = [i for i in range(b)]
    A.remove(d)
    N = b - 1

    lowblocks = [[0]]  # filler

    A1 = [a for a in A if a !=0]

    lowblocks.append(A1)
    # print(A1)

    maxblock = []
    maxblock.append([b * x + a for x in A1 for a in A])
    if d == 0:
        maxblock.append([b * i for i in range(1, b)])
        maxblock.append([])
    else:
        L = [b * a + d for a in A1]
        L.extend([b * d + a for a in A])
        maxblock.append(L)
        maxblock.append([b * d + d])

    if verbose:
        print("cardinalité au niveau 2 pour k=0, 1, 2 (avec chiffre %s): %s, %s, %s" %
              (d, len(maxblock[0]), len(maxblock[1]), len(maxblock[2])))

    assert len(maxblock[0]) + len(maxblock[1]) + len(maxblock[2]) == b*(b-1), "Fix this block2 bug!"

    lesgammas = [ N ]
    for j in range(1, Mmax+1):
        lesgammas.append(float(sum(a**j for a in A1)))
    
    lespuissancesded = [ 1 ]
    for j in range(1, Mmax+1):
        lespuissancesded.append(float(d**j))

    # calcul récursif des moments
    touslescoeffs = []

    # k = 0
    coeffs = [ float(b) ]
    for m in range(1, Mmax+1):
        coeffs.append(sum(Pascal[m][j]*lesgammas[j]*coeffs[m-j] for j in range(1,m+1))
                      / (b**(m+1)-N)  # N = b-1
                      )  

    touslescoeffs.append(coeffs)

    for j in range(1,k+1):
        prevcoeffs = coeffs

        coeffs = [ float(b) ]
        for m in range(1, Mmax+1):
            coeffs.append((sum(Pascal[m][j]*lesgammas[j]*coeffs[m-j] for j in range(1,m+1))
                           + sum(Pascal[m][j]*lespuissancesded[j]*prevcoeffs[m-j] for j in range(0,m+1))
                           )
                          / (b**(m+1)-N))  # N = b-1
        touslescoeffs.append(coeffs)
   
    # boucle pour évaluer si all = True également les j < k

    for j in range(0 if all else k, k+1):

        # calcul de la série alternée de Burnol

        S = 0

        if j == 0:
            S = sum(1/x for x in A1)

        if j == 1:
            if d != 0:
                S = 1/d

        if verbose:
            print("sommes du niveau 1 pour d = %s et j = %s: %s" % (d, j, S))

        if j >= 2:
            S += b * sum(1/x for x in range(b, b**2))
        else:
            S += b * sum(1/x for x in maxblock[0])
            if j == 1:
                S += b * sum(1/x for x in maxblock[1])

        if verbose:
            print("somme ajustée de niveau 2 avant corrections géométriques: %s" % S)

        #  print(Mmax, b, len(maxblock), N, l, f)
        if verbose:
            print("on va utiliser %d termes de la série alternée" % Mmax)

        bubu = touslescoeffs[j][Mmax] * sum(1/(n ** (Mmax+1)) for n in maxblock[0])
        if j >= 1: 
            bubu += touslescoeffs[j-1][Mmax] * sum(1/(n ** (Mmax+1)) for n in maxblock[1])
        if j >= 2:
            bubu += touslescoeffs[j-2][Mmax] * sum(1/(n ** (Mmax+1)) for n in maxblock[2])

        if verbose:
            print("le plus petit terme retenu est %s" % (-bubu if (Mmax & 1) else bubu))

        for m in range(Mmax-1, 0, -1):  # last one is m=1
            bubu = -bubu
            bubu += touslescoeffs[j][m] * sum(1/(n ** (m+1)) for n in maxblock[0])
            if j >= 1:
                bubu += touslescoeffs[j-1][m] * sum(1/(n ** (m+1)) for n in maxblock[1])
            if j >= 2:
                bubu += touslescoeffs[j-2][m] * sum(1/(n ** (m+1)) for n in maxblock[2])

        if verbose:
            print("le terme correctif qui sera retiré est: %s" % bubu)
            print("b = %s, d = %s, j = %s" % (b, d, j))

        S = S - bubu

        if all:
            print("(k=%2d) %s" % (j, S))

    return S


def irwinposfloat(b, d, k, all=False, verbose=False, Mmax = -1):
    """Calcule la somme d'Irwin pour la base b et le chiffre d et l'entier k

    Utilise algorithme de Burnol, série positive de niveau 2

    :param b: int
        la base
    :param d: int
        le chiffre
    :param k: int
        l'entier k
    :param all: bool,optional (False)
        si True, on calcule et on affiche
        également pour tous les entiers inférieurs ou égaux à k.
    :param verbose: bool,optional (False)
        si True affiche informations supplémentaires
    :param Mmax: int,optional (-1)
        si positif ou nul, donne le nombre de termes de la série de Burnol
        à utiliser
    :rtype: float
    :return: la somme d'Irwin de hauteur k pour le chiffre d en base b

    """

    assert b>1, "%s doit être au moins 2" % b

    assert 0 <= d < b, "%d doit être positif et au plus b-1" % d

    if Mmax == -1:
        Mmax = ceil(15/log10(b)) +2  # +2 de sécurité utile que pour des grands b
                                     # mais pour ceux-là de toute façon algo pas praticable
                                     # le +2 est clairement overkill
    if verbose:
        print("Mmax = ", Mmax)

    # # auxiliaire
    # def fb(a):
    #     return lambda x : b * x + a

    A = [i for i in range(b)]
    A.remove(d)
    # print(A)
    N = b - 1

    lowblocks = [[0]]  # filler

    A1 = [a for a in A if a !=0]

    lowblocks.append(A1)
    # print(A1)

    maxblock = []
    maxblock.append([b * x + a for x in A1 for a in A])
    if d == 0:
        maxblock.append([b * i for i in range(1, b)])
        maxblock.append([])
    else:
        L = [b * a + d for a in A1]
        L.extend([b * d + a for a in A])
        maxblock.append(L)
        maxblock.append([b * d + d])

    if verbose:
        print("cardinalité au niveau 2 pour k=0, 1, 2 (avec chiffre %s): %s, %s, %s" %
              (d, len(maxblock[0]), len(maxblock[1]), len(maxblock[2])))

    assert len(maxblock[0]) + len(maxblock[1]) + len(maxblock[2]) == b*(b-1), "Fix this block2 bug!"

    # ATTENTION que la série positive a des récurrences avec b-1-d à la place de d
    Aprime = [i for i in range(b)]
    dprime = b-1-d
    Aprime.remove(dprime)

    A1prime = [a for a in Aprime if a != 0]
    lesgammasprime = [ N ]
    for j in range(1, Mmax+1):
        lesgammasprime.append(float(sum(a**j for a in A1prime)))
    
    lespuissancesdedprime = [ 1 ]
    for j in range(1, Mmax+1):
        lespuissancesdedprime.append(float(dprime**j))

    # calcul récursif des moments
    touslescoeffs = []

    # k = 0
    coeffs = [ float(b) ]
    for m in range(1, Mmax+1):
        coeffs.append((b**(m+1)  # terme additionnel pour la récurrence des v
                       + sum(Pascal[m][j]*lesgammasprime[j]*coeffs[m-j] for j in range(1,m+1))
                       )/ (b**(m+1)-N)  # N = b-1
                      )  

    touslescoeffs.append(coeffs)

    for j in range(1,k+1):
        prevcoeffs = coeffs

        for m in range(1, Mmax+1):
            coeffs = [ float(b) ]
            for m in range(1, Mmax+1):
                coeffs.append((sum(Pascal[m][j]*lesgammasprime[j]*coeffs[m-j] for j in range(1,m+1))
                               + sum(Pascal[m][j]*lespuissancesdedprime[j]*prevcoeffs[m-j] for j in range(0,m+1))
                               )
                              / (b**(m+1)-N))  # N = b-1
        touslescoeffs.append(coeffs)
   
    # boucle pour évaluer si all = True également les j < k

    for j in range(0 if all else k, k+1):

        # calcul de la série poistive de Burnol

        S = 0

        if j == 0:
            S = sum(1/x for x in A1)

        if j == 1:
            if d != 0:
                S = 1/d

        if verbose:
            print("sommes du niveau 1 pour d = %s et j = %s: %s" % (d, j, S))

        if j >= 2:
            S += b * sum(1/(x + 1) for x in range(b, b**2))
        else:
            S += b * sum(1/(x + 1) for x in maxblock[0])
            if j == 1:
                S += b * sum(1/(x + 1) for x in maxblock[1])

        if verbose:
            print("somme ajustée de niveau 2 avant corrections géométriques: %s" % S)

        #  print(Mmax, b, len(maxblock), N, l, f)
        if verbose:
            print("on va utiliser %d termes de la série positive" % Mmax)

        bubu = touslescoeffs[j][Mmax] * sum(1/((n + 1) ** (Mmax+1)) for n in maxblock[0])
        if j >= 1: 
            bubu += touslescoeffs[j-1][Mmax] * sum(1/((n + 1) ** (Mmax+1)) for n in maxblock[1])
        if j >= 2:
            bubu += touslescoeffs[j-2][Mmax] * sum(1/((n + 1) ** (Mmax+1)) for n in maxblock[2])

        if verbose:
            print("le plus petit terme retenu est %s" % bubu)

        for m in range(Mmax-1, 0, -1):  # last one is m=1
            bubu += touslescoeffs[j][m] * sum(1/((n + 1) ** (m+1)) for n in maxblock[0])
            if j >= 1:
                bubu += touslescoeffs[j-1][m] * sum(1/((n + 1) ** (m+1)) for n in maxblock[1])
            if j >= 2:
                bubu += touslescoeffs[j-2][m] * sum(1/((n + 1) ** (m+1)) for n in maxblock[2])

        if verbose:
            print("le terme correctif qui sera ajouté est: %s" % bubu)
            print("b = %s, d = %s, j = %s" % (b, d, j, level))

        S = S + bubu

        if all:
            print("(k=%2d) %s" % (j, S))

    return S

def irwinfloattest():
    """Comparaison avec données qu'on peut récupérer dans article de Baillie, 2008
    """
    # données de Baillie, extraites de son arxiv de 2008 car je n'ai pas
    # de Mathematica donc je ne peux pas faire tourner son programme pour
    # comparer et ne peux utiliser que les données qu'on trouve dans
    # les exemples listés

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

    # somme avec 3  zéros : 23.025851037148538
    # somme avec 10 zéros : 23.025850929940457
    # somme avec 0 zéro en base 2: 1.60669 51524 15291 76378
    # somme avec 1 zéro en base 2: 1.46259 07350 44364 69954 61454

    for d in range(10):
        for k in range(3):
            x = irwinfloat(10,d,k); print(x, " (%s, %s, %s)" % (10, d, k));
            delta = x - bailliedata[d][k]; print(delta)
            assert abs(delta) < 2e-14,\
                "problème avec b=10, d=%s, k=%s, série alternée" % (d, k)

            x = irwinposfloat(10,d,k); print(x, " (%s, %s, %s) (pos)" % (10, d, k));
            delta = x - bailliedata[d][k]; print(delta)
            # la série converge moins bien que la série alternée
            assert abs(delta) < 3e-14,\
                "problème avec b=10, d=%s, k=%s, série positive" % (d, k)

    for data in ([(10,0,3), 23.025851037148538],
                 [(10,0,10),23.025850929940457],
                 [(2,0,0),   1.60669515241529176378],
                 [(2,0,1),   1.46259073504436469954],
                 [(2,1,0),   0.],  #  faudrait tester égalité mais bon de toute façon ok
                 [(2,1,1),   2.]): #  pas la peine de compliquer inutilement
        assert abs(irwinfloat(*data[0]) - data[1]) < 2e-14,\
            "problème avec b=%S, d=%s, k=%s, série alternée" % (data[0][0],
                                                                data[0][1],
                                                                data[0][2])
        assert abs(irwinposfloat(*data[0]) - data[1]) < 3e-14,\
            "problème avec b=%S, d=%s, k=%s, série positive" % (data[0][0],
                                                                data[0][1],
                                                                data[0][2])


if __name__ == "__main__":
    print("""Utilisez l'option -i de python si python quitte immédiatement!
    
irwinfloat(b, d, k):
   calcule suivant l'algorithme de Burnol la somme d'Irwin des
   1/n pour n ayant en base "b" le chiffre "d" présent exactement "k" fois.

   La valeur est calculée avec des floats et les erreurs d'arrondis rendent
   au moins le dernier chiffre sans signification.

irwinfloat(b, d, k, all = True)
   calcule et affiche non seulement pour k mais pour tous les j entre 0
   et k

irwinposfloat(b, d, k) et irwinposfloat(b, d, k, all=True) sont des variantes
qui utilisent la série positive de Burnol et non pas sa série alternée.

*Dans les deux cas, le dernier chiffre et même les deux derniers sont sans
signification.*"""
    )
