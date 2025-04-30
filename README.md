# irwin


## Description

This repository is to document implementations of the formulas from the
February 2024 article
[Measures for the summation of Irwin series](https://arxiv.org/abs/2402.09083)
by me (aka Jean-François Burnol).

## Files

So far all contributed files with extension `.sage` are to be used
interactively within a SageMath session.  Start a SageMath session and at the
`sage` prompt, enter `load("<FOO>.sage")`.  Follow banner instructions which
will be printed.  For example:

```sage
sage: load("irwin_v3.sage")
[some info printed]
sage: irwin(10,9,0,52)
22.92067661926415034816365709437593191494476243699848
```

Hundreds, even thousands of digits, can be computed.  For this, prefer using
`irwin_v5.sage` which will try to take advantage of multiple cores on your
system.

```sage
sage: load("irwin_v5.sage")
[some info printed mentioning maxworkers which defaults to 8]
sage: irwin(10,9,0,1002)
[the output is reformatted here to show 50 decimals per line]
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
```

We can also compute the Irwin harmonic sums, where a given digit is allowed
only a given number of time:

```sage
sage: irwin(10,9,1,1002)
[the output is reformatted here to show 50 decimals per line]
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
```

Thanks to Arnaud Bodin and Yusuf Emin Akpınar for sharing their interest into
the Irwin paper and its predecessor
[Moments in the exact summation of the curious series of Kempner type](https://arxiv.org/abs/2402.08525),
and in particular for reporting progress in computing tens of thousands of
digits (Y. E. A. has obtained 100000=2+99998 digits of the classical "no-9"
radix 10 Kempner series).  This motivated me to revisit the 2024 code with the
perspective to make it easier and less demanding on the hardware to obtain
10000+ digits even within some high-level front-end such as SageMath.

- [irwinfloat.py](irwinfloat.py) and [irwin.sage](irwin.sage) are the files
  available at [my arXiv paper](https://arxiv.org/abs/2402.09083). They were
  last updated there at v4 of April 6, 2024.

- [irwin_v3.sage](irwin_v3.sage) is an evolution ot the latter which is better
  suited to computing thousands of digits (say, starting at around 2000
  digits).  It considerably reduces the memory footprint, and has some speed
  gain from using multiple precisions for the `RealField` objects.  On a
  current personal computer it will obtain 50000 digits of the classic Kempner
  sum in a matter of hours (depending on your hardware).  This computation
  attempted with the 2024 version would have needed (with level=3, hence
  needing binomial coefficients up to `N=25000`) **at least 470 gigabytes** of
  computer memory.  In contrast, the 2025 version needs (theoretically,
  i.e. using the base-2 logarithm as estimate of the number of bits needed) of
  the order of 113 megabytes.

  The `irwin()` and `irwinpos()` provided by this file use by default
  `level=3` which is better for getting thousands of digits.  They also have a
  slightly safer way of estimating how many terms of the series should be kept
  in order for the result to be the correctly rounded one, also if using
  `level=4` or if handling Irwin series (`k>0`). The code does not check
  though if the computed value is near to be exactly half-way between floor
  and ceil, so to be absolutely certain the printed value is really the
  correctly rounded one, one should redo the computation asking for some extra
  digits (say 5 more)...

- [irwin_v3_loader.py](irwin_v3_loader.py) and
  [irwin_loader.py](irwin_loader.py) and similarly named files can be used to
  load concurrently [irwin.sage](irwin.sage) and
  [irwin_v3.sage](irwin_v3.sage) into the same interactive SageMath session.
  See the comments therein.  No guarantees though.

- [irwin_v5.sage](irwin_v5.sage) uses "parallelization".  This is very
  efficient for the computation of the `beta(m+1)`'s.  It is more difficult to
  apply parallelization to the computation of the recurrent coefficients
  `u_{k;m}` (or `v_{k;m}` for the series with positive coefficients).  In an
  earlier `v4` version, now obsolete, this was done coefficient per
  coefficient, via a division of its defining sum into subsums.  The `v5`
  version proves more efficient: it does not use subsums but computes in
  parallel for successive indices `m=M+1`, ..., `m=M+maxworkers`, up to
  finitely many additive corrections done after the parallelized subprocesses
  have returned; this is done for all `j`'s from `0` to `k` in one-go.

  The memory foot-print is higher than in `v3` because `maxworkers` rows of
  the Pascal triangle are in memory rather than only `2` for `v3` (and `v4`).

  The variable `maxworkers` does not have to be at most the actual number of
  cores on the user system.  We tested both `maxworkers=8` and `maxworkers=2`
  on an old hardware with only 2 cores, timings were about the same.

  **Keep in mind any change to `maxworkers` must be followed by a re-load of
  `irwin_v5.sage` in the interactive session.**

- I should note here that there are some surprises with parallelization on my
  macOS 15.4.1 Sequoia.  To parallelize the computation of the moments of the
  measure defined in my research, for which the known way is to proceed via
  the recurrences I obtained, the code necessarily has to call many times
  procedures which have been `@parallel` decorated; for example with
  `maxworkers=8` and if we need `10000` such coefficients, we will call `1250`
  times such a procedure.  Turns out that on macOS 15.4.1 this causes a linear
  increase in some sort of waiting time which increases the execution times.
  It is hard to analyse on `irwin_v5.sage` because the time increases at least
  linearly with the index `m` purely due to the length of the recurrence (as
  one sees with the test timings reported regularly if with the
  `showtimes=True` option).  Anyway I reduced this to a reproducer
  [test_parallel_sleep.sage](test_parallel_sleep.sage) and I reported as
  [this SageMath ticket](https://github.com/sagemath/sage/issues/39960) which
  has been closed, probably more as a "wont-fix" (see also #1).  I consider
  the ticket was closed (twice actually) a bit fast, with no obvious desire to
  investigate more what appears to be specific to macOS 15.4.1.  I have some
  difficulties to fathom how a more powerful processor than most others could
  behave so bad, and it definitely is cause to ponder why the caching done by
  the OS has such an obvious detrimental effect in this context.  Is it a bug
  somewhere in the OS or is it a bug in Python?  It is hard for me to consider
  it a fact of life that when you call many times such a procedure it affects
  what is done next... unless you kill `sage` and reload! I must admit I am
  saying this under the impression that the problem does not show on Linux
  boxes but I did not investigate fully myself for lack of resources.  All I
  know is that on a quadri-core much older macOS there is absolutely **no**
  problem.  Anyway, I am quite interested into reports on how
  [test_parallel_sleep.sage](test_parallel_sleep.sage) behaves on various
  systems.  Please comment at #1.

  In practice this means that on macOS 15 at least one should quit `sage` and
  restart a session before launching a run for which one needs the best
  efficiency.  The `irwin_v5` module (version 1.5.4 of April 23, 2025) has
  computed via `irwin()` on a Mac Mini M4 Pro with `10+4` cores `2+101010`
  decimals of the "no-9" Kempner constant in 5h38mn, with `maxworkers` left to
  its default `8`.
  
- The next files with names of the type `k_prec_N` contain decimal expansions
  of the classic "no-9 radix-10" Kempner series `22.92067661926415...`,
  correctly rounded to `N` decimal places.  They are currently:
  * [k_prec_1000](k_prec_1000)
  * [k_prec_2000](k_prec_2000)
  * [k_prec_5000](k_prec_5000)
  * [k_prec_10000](k_prec_10000)
  * [k_prec_20000](k_prec_20000)
  * [k_prec_25000](k_prec_25000)
  * [k_prec_50000](k_prec_50000)
  * [k_prec_99998](k_prec_99998) (contributed by Yusuf Emin Akpınar, thanks!)
  
- [taille_pascal.pdf](taille_pascal.pdf) explains how many bits are needed to
  store in computer memory the Pascal triangle up (or rather down) to a
  certain row, or the memory needed to store only one such row.  Thanks to
  Nicolas Radulesco who asked for SageMath equivalents of the few Maple usages
  which were inserted in the text, they are in
  [taille_pascal_symbolic](taille_pascal_symbolic).

## TODO

Add some additional bibliographical references.

## License

The files in this repository are distributed under the
CC-BY-SA 4.0 License.  See [LICENSE](LiCENSE).
