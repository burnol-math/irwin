================
burnolmath/irwin
================

This page serves as a mathematical introduction to
https://gitlab.com/burnolmath/irwin.  This latter repository aim is to provide
the software and to comment on implementation aspects, but not to go into
mathematics except when this is relevant to the implementation.  Our account
will be brief.  Follow the links to the original literature for more details.

.. contents::

The Kempner and Irwin harmonic series
=====================================

Let $b>1$ be an integer (it will be called *base* or *radix*),
$d\in\{0,...,b-1\}$ a digit, and $k$ a non-negative integer.  We associate to
them $S(b,d,k) = \sum\frac{1}{n}$ where only those positive integers $n$ are
kept whose radix-$b$ (aka $b$-ary) representation contains *exactly* $k$
occurrences of the digit $d$.

\F. Irwin established the convergence of this sub-series in 1916 in `A Curious
Convergent Series <irwin_>`_.  Earlier in 1914, A. J. Kempner had initiated
the topic in a note with the same title `A Curious Convergent Series
<kempner_>`_.  He made the observation that keeping from the harmonic series
only those terms not using the decimal digit $9$ (or any other given decimal
digit) gives a convergent series.  Both Kempner and Irwin considered only
$b=10$.  In the special case with $b=2$, $d=1$ and $k=0$, $S(2,1,0)$ is
defined to be zero, as it is the value of an empty sum.

.. _irwin: https://doi.org/10.2307/2974352
.. _kempner: https://doi.org/10.2307/2972074

A simple proof of convergence in the "no-$9$" Kempner case goes as follows:
among integers $n$ verifying $10^{l-1}\leq n< 10^l$, i.e. requiring $l$ digits
exactly, the first one not zero, there are exactly $8\cdot 9^{l-1}$ having no
occurrence of $9$ among their decimal digits.  The combined contribution of
all those $\frac{1}{n}$'s to the harmonic series is thus bounded above by
$8\cdot 9^{l-1} \cdot 10^{-(l-1)}$ and $S(10,9,0)<\sum_{l=1}^\infty 8\cdot
(0.9)^{l-1} = 8/(1 - 0.9) = 80$. The same reasoning applies to any other
excluded digit (but gives $90$ as upper bound if the excluded digit is $0$).

We can improve the upper bound of the previous paragraph easily to
$\sum_{d=1}^8\frac{1}{d} + 8\frac{0.9}{1 - 0.9} < 75$.  But this is still far
from the actual value.  Via a more clever argument Irwin showed
$S(10,9,0)<23.3$ and also established $22.4$ as lower bound.  One could think
that obtaining lower bounds is an easy task as any partial sum of the series
will serve at such.  This is not practical though because the convergence is
very slow: as explained by R. Baillie (link below) the sum of all terms in the
series with denominators having at most 30 digits is still less than $22$.  So
Irwin did not obtain the $22.4$ lower bound simply by summing enough many
terms...

This page is not meant to be encyclopedic on the topic so we do not give here
the details of Irwin reasoning.

Irwin actually considered more generally the following condition: only keep
those $n$'s with at most a given number of $0$'s, *or* at most a given number
of $1$'s, *or* at most a given number of $2$'s, ..., and he proved the
convergence in that case.  The convergence reduces easily to the convergence
of the $S(b,d,k)$'s.  We limit our consideration to the latter here because they
are the series for which the repository https://gitlab.com/burnolmath/irwin
currently provides numerical implementations.


The Baillie algorithm
=====================

In 1979, Baillie proposed a method to compute accurate numerical values of the
Kempner sums $S(10,d,0)$ for $d=0$, $1$, ..., $9$.  Here they are, truncated
to 20 decimal places as they were published in the Baillie paper::

    S(10,0,0) = 23.10344 79094 20541 61603...
    S(10,1,0) = 16.17696 95281 23444 26657...
    S(10,2,0) = 19.25735 65328 08072 22453...
    S(10,3,0) = 20.56987 79509 61230 37107...
    S(10,4,0) = 21.32746 57995 90036 68663...
    S(10,5,0) = 21.83460 08122 96918 16340...
    S(10,6,0) = 22.20559 81595 56091 88416...
    S(10,7,0) = 22.49347 53117 05945 39817...
    S(10,8,0) = 22.72636 54026 79370 60283...
    S(10,9,0) = 22.92067 66192 64150 34816...

We do not describe here the details of Baillie's method, which are explained
in the following references:

- Robert Baillie: Sums of reciprocals of integers missing a given
  digit. Amer. Math. Monthly 86(5), 372–374 (1979). `DOI link`__

  __ https://doi.org/10.2307/2321096
- Robert Baillie: Summing the curious series of Kempner and Irwin (2008).
  `arXiv:0806.4410 <https://arxiv.org/abs/0806.4410>`_.
- Thomas Schmelzer and Robert Baillie: Summing a curious, slowly convergent
  series. Amer. Math. Monthly 115(6), 525–540 (2008). `DOI link`__

  __ https://doi.org/10.1080/00029890.2008.11920559

In the latter paper of Schmelzer and Baillie an interesting further
generalization is considered about excluding not a digit but a *word* of
multiple contiguous digits.  This will not be touched upon here.  But let us
mention that Allouche, Hu and Morin have shown the interesting result that for
such a word $w$ of $p$ digits in base $b$, and a number of occurrences $k$
going to infinity, the corresponding generalized Irwin harmonic series
converge to $b^p\log(b)$ (with the natural logarithm here of course).  For
$p=1$ this had been shown earlier by Farhi.

- Farhi, B., A curious result related to Kempner’s series. Amer. Math. Monthly
  115(10), 933–938 (2008). `DOI link`__

  __ https://doi.org/10.1080/00029890.2008.11920611
- Allouche, J.-P., Hu, Y. and Morin, C.: Ellipsephic harmonic series
  revisited.  Acta Math. Hungar. 173, 461–470 (2024). `DOI link`__

  __ https://doi.org/10.1007/s10474-024-01448-5

The new series
==============

In 2024, the author of this page published novel theoretical formulas
representing exactly the Kempner-Irwin series.  Here are these formulas, as
proven in `Measures for the summation of Irwin series
<https://arxiv.org/abs/2402.09083>`_ and earlier for $k=0$ (and only the
alternating series) in `Moments in the exact summation of the curious series
of Kempner type <https://arxiv.org/abs/2402.08525>`_.

Let $b>1$, $d$ a $b$-ary digit, and $k$ a non-negative integer.  We will need
the following notation: for $n$ a positive integer, $k(n)$ is the number of
occurrences of digit $d$ in the $b$-ary representation of $n$.  The notation
drops mention of $b$ and $d$ to avoid becoming too heavy.

Define first these power sums over the digits other than $d$ (here too we drop
$b$ and $d$ from the notation):

$$\gamma_n = \sum_{a\neq d,\; 0\leq a<b} a^n$$

Then define (rational) coefficients $u_{j;m}$ via the following recurrences:

- $u_{0;0} = b$ and for $m\geq1$:

  $$u_{0;m} = (b^{m+1} - b +1)^{-1}\sum_{n=1}^{m} \binom{m}{n} \gamma_n u_{0:m-n}$$
- For $j\geq1$, $u_{j;0}=b$ and for $m\geq1$:
  $$u_{j;m}=(b^{m+1} - b +1)^{-1}\left(\sum_{n=1}^{m}\binom{m}{n}\gamma_n u_{j;m-n}
  + \sum_{n=0}^{m}\binom{m}{n}d^n u_{j-1;m-n}\right)$$

Let now $l$ be a positive integer (also called the "level").  Preferably we will
take $l$ at least $2$ to guarantee that the series defined next has geometric
convergence.

We then have the following exact formula, which uses the $k(n)$ notation
defined earlier:

$$S(b,d,k) = \sum_{1\leq n<b^{l-1}, k(n)=k}\frac1n
+ b\sum_{b^{l-1}\leq n<b^l, k(n)\leq k}\frac1n
+ \sum_{m=1}^\infty (-1)^{m}\sum_{b^{l-1}\leq n<b^l, k(n)\leq k}
\frac{u_{k-k(n);m}}{n^{m+1}}$$

It can be shown that the alternating series has its contributions decreasing
in absolute values (except if $b=2$, $d=1$, $k=0$ or $k=1$ as then all
contributions to the alternating series vanish).  So using its partial sums
gives upper and lower bounds.  For example dropping the alternating series
altogether (whose first term is negative) gives as upper bound for
$S(10,9,0)$:

$$S(10,9,0) < \sum_{1\leq d\leq 8}\frac1d + 10\sum_{1\leq d\leq 8, 0\leq e\leq
8}\frac1{10d+e} =23.2577...$$

(the numerical value is truncated). It is interesting that this is actually
exactly the upper bound given by Irwin in his paper, but he computed
it numerically with only one decimal place of precision and obtained $23.3$ as
upper bound.

There is another exact formula, which uses only positive terms.
Define first these power sums over the digits other than $b-1-d$:

$$\gamma_n' = \sum_{a\neq b-1-d,\; 0\leq a<b} a^n$$

Then define (rational) coefficients $v_{j;m}$ via the following recurrences:

- $v_{0;0} = b$ and for $m\geq1$:

  $$v_{0;m} = (b^{m+1} - b +1)^{-1}\left(b^{m+1}
  +\sum_{n=1}^{m} \binom{m}{n} \gamma_n' v_{0:m-n}\right)$$
- For $j\geq1$, $v_{j;0}=b$ and for $m\geq1$:
  $$v_{j;m}=(b^{m+1} - b +1)^{-1}\left(\sum_{n=1}^{m}\binom{m}{n}\gamma_n' v_{j;m-n}
  + \sum_{n=0}^{m}\binom{m}{n}(b-1-d)^n v_{j-1;m-n}\right)$$

Let now $l$ be positive integer (also called the "level").
We then have this further exact formula:

$$S(b,d,k) = \sum_{1\leq n<b^{l-1}, k(n)=k}\frac1n
+ b\sum_{b^{l-1}\leq n<b^l, k(n)\leq k}\frac1{n+1}
+ \sum_{m=1}^\infty \sum_{b^{l-1}\leq n<b^l, k(n)\leq k}
\frac{v_{k-k(n);m}}{(n+1)^{m+1}}$$

Dropping the series indexed by $m\geq1$ altogether (all whose terms are
positive, except for $b=2$, $d=1$, $k=0$, then they all vanish) gives the
following lower bound for $S(10,9,0)$:

$$S(10,9,0) > \sum_{1\leq d\leq 8}\frac1d + 10\sum_{1\leq d\leq 8, 0\leq e\leq
8}\frac1{10d+e+1}= 22.4249\dots$$

It is interesting that Irwin had actually obtained a slightly sharper lower bound:

$$S(10,9,0) > \sum_{1\leq d\leq 8}\frac1d + \sum_{1\leq d\leq 8, 0\leq e\leq
8}\frac1{10d+e} + 9\sum_{1\leq d\leq 8, 0\leq e\leq 8}\frac1{10d+e+1}
=22.5081\dots $$

But he computed it with not enough precision and ended the paper with the
lower bound $22.4$ where he could have claimed $22.5$.  This Irwin lower
bound can be deduced from the series using the fact that in this case
the $v_{0;m}$'s are larger than $1$ (they are even $>10/9$).


SageMath implementation
=======================

The formulas of the previous section have been implemented in a SageMath_
script ``irwin.sage`` available at https://gitlab.com/burnolmath/irwin.

.. _SageMath: https://www.sagemath.org

Here is an example of use::

    sage: load("irwin.sage")
    [some info printed]
    sage: for d in range(10):
    ....:     print(f"S(10,{d},0) = {irwin(10,d,0,52)}")
    ....: 
    S(10,0,0) = 23.10344790942054161603405404332559813830280000528214
    S(10,1,0) = 16.17696952812344426657960388036400930556721979076313
    S(10,2,0) = 19.25735653280807222453277677019445411552605383115487
    S(10,3,0) = 20.56987795096123037107521741905311141415386967473078
    S(10,4,0) = 21.32746579959003668663940148693951284375095170327002
    S(10,5,0) = 21.83460081229691816340723504060918271784656751501392
    S(10,6,0) = 22.20559815955609188416738048000752710519385610666846
    S(10,7,0) = 22.49347531170594539817622691533977597400591554167251
    S(10,8,0) = 22.72636540267937060283364415674255788921070261636022
    S(10,9,0) = 22.92067661926415034816365709437593191494476243699848
    sage: for d in range(10):
    ....:     print(f"S(10,{d},1) = {irwin(10,d,1,52)}")
    ....: 
    S(10,0,1) = 23.02673534156912696109462698601416425917373603671403
    S(10,1,1) = 23.16401859427283204084669788222982114096553545927615
    S(10,2,1) = 23.08826066275634239334138412149393536360742700309969
    S(10,3,1) = 23.06741088193023010241930360133553872732553405199013
    S(10,4,1) = 23.05799241338182439575664489721510980557994895159896
    S(10,5,1) = 23.05272889453011749903870693765255517445588205770562
    S(10,6,1) = 23.04940997329550055703597736374339006623306171815468
    S(10,7,1) = 23.04714619019864185082971899384931544629351115531788
    S(10,8,1) = 23.04551390798215553341865446180006252186687724092286
    S(10,9,1) = 23.04428708074784831967594930973617482538959203064774

Hundreds or thousands of decimal digits are computed easily::

    sage: load("irwin.sage")
    [some info printed mentioning maxworkers which defaults to 8]
    sage: irwin(10,9,0,1002)
    "the output is reformatted here to show 50 decimals per line"
    22.92067661926415034816365709437593191494476243699848
       15685419983565721563381899111294456260374482018989
       90964125332346922160471190478310297506146968857121
       01806786493339402886962779578685961198637905620169
       32188040880170136179021106286611735099211021080576
       70378581471208344258765832272657620103831470760370
       30815999623544735896526905676888497081960327431233
       14588927997290413878495249814944204592152773507367
       07218520004083026308916169121123862636859589823575
       17170592498667879488473210892480659162340101523560
       00506548043749678309013031335561096953014813317749
       55762523805629716085009843545476018253422157510734
       48392165782984461954239160106117835383539414385364
       56085452218993239443664387904158857609144227813991
       99992242055353569500690341681751890944480911928277
       83446999651712608600666360667788028808406885936480
       28751790909188136795127797348003365941380076337136
       20275923523021897838806069615932191066192832138116
       95786715012908593756769518010810881852946961772722
       23692633510303284693132263332046629826719621921950
    sage: irwin(10,9,1,1002)
    "the output is reformatted here to show 50 decimals per line"
    23.04428708074784831967594930973617482538959203064773
       62135578783008262042579280261007145671482118830782
       57921943364250671219896744152423420975036205340023
       04185784253363785559895004750435973037013624248692
       81693758690484874474533229028198473110373901313573
       11432258255340503360367871022252892752126397786453
       83974571210767152137049460532051975462039681885418
       36330606303719288329589891730925547573307176245910
       08049079306964817292360040758680171370594059804295
       02478333846389429682664501235076520170203860074879
       38353118391668568399204076846814722810630157329244
       15926574398873339404075981049839733114335471459356
       69038199984752130884656309345224669665996152072319
       08005559300473069297421785342383181330939055266000
       60305950156550830604733577855412791540396401495913
       50606434552699589093508098625808576790394619236238
       07758802238788028673496735651848096113580106792051
       15685105510116478927387735732377302778507573807355
       85556733655125237297349752705490397573629470971925
       89207555878085304702890838585263137514388675244390


Use ``irwinpos()`` for the implementation of the positive series rather than
the alternating one.
Check ``help(irwin)`` or ``help(irwinpos)`` for additional parameters.


References
==========

For $k=0$ the alternating series was first published in:

- Moments in the exact summation of the curious series of Kempner type,
  https://arxiv.org/abs/2402.08525.

The latter paper (which is to appear in *Amer. Math. Monthly* in 2025 or 2026)
also gives exact theoretical formulas for Kempner-like series where multiple
digits are simultaneously excluded.

The alternating and positive series for all $k\geq0$ have been first
published in:

- Measures for the summation of Irwin series, https://arxiv.org/abs/2402.09083.

This manuscript is currently submitted.  Some of my further research
has already appeared:

- Summing the "exactly one 42" and similar subsums of the harmonic series,
  Advances in Applied Mathematics Volume 162, January 2025, 102791. `DOI link`__

  __ https://doi.org/10.1016/j.aam.2024.102791

- Digamma function and general Fischer series in the theory of Kempner sums,
  Expositiones Mathematicae, Volume 42, Issue 6, December
  2024, 125604. `DOI link`__

  __ https://doi.org/10.1016/j.exmath.2024.125604

- Measures associated with certain ellipsephic harmonic series and the
  Allouche-Hu-Morin limit theorem, Acta Mathematica Hungarica (2025).
  `DOI link`__

  __ https://doi.org/10.1007/s10474-025-01525-3

These next two manuscripts by the author are awaiting referee reports.

- Sur l'asymptotique des sommes de Kempner pour de grandes bases,
  `arXiv:2403.01957 <https://arxiv.org/abs/2403.01957>`_.
- Un développement asymptotique des sommes harmoniques de Kempner-Irwin,
  `arXiv:2404.13763 <https://arxiv.org/abs/2404.13763>`_.
