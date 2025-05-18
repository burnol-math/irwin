# irwin_loader.py
# ceci permet de faire depuis l'invite de sage, par exemple:

# from irwin_loader import irwin as irwin5, irwinpos as irwinpos5

# Ceci par exemple
# afin de comparer avec le comportement des procédures
# du irwin_legacy.sage de 2024 ou de irwin_v3.sage de 2025.

# Cette méthode a de très gros problèmes potentiels
# je pense si des fonctions auxiliaires de mêmes noms
# mais différentes sont dans les fichiers.
# Le souci est que je n'ai jamais pu accéder à
# sage.misc.preparser.  Mais je pense que pour
# comparer irwin_v3.sage (2025) et irwin.sage (2024)
# ça va marcher.

from sage.all import *

load("irwin.sage")
irwin = globals()['irwin']
irwinpos = globals()['irwinpos']
