# irwin_loader.py
# ceci permet de faire depuis l'invite de sage, par exemple:

# from irwin_loader import irwin as irwin2024, irwinpos as irwinpos2024

# Ceci
# afin de comparer avec le comportement
# du irwin.sage de 2024 qui définit des fonctions
# du même nom.

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
