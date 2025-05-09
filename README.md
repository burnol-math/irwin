# irwin

[[_TOC_]]

## Description

This repository is to document implementations of the formulas from the
February 2024 article
[Measures for the summation of Irwin series](https://arxiv.org/abs/2402.09083)
by me (aka Jean-François Burnol).

## Quick guide

The main files currently implementing the algorithm are `irwin_v3.sage` and
`irwin_v5.sage`.  Start a SageMath session and at the `sage` prompt, enter
`load("irwin_v3.sage")` or `load("irwin_v5.sage")`.  Follow banner
instructions which will be printed.  For example:

```sage
sage: load("irwin_v3.sage")
[some info printed]
sage: irwin(10,9,0,52)
22.92067661926415034816365709437593191494476243699848
```

For computing Kempner harmonic sums with thousands of decimal digits, or Irwin
sums starting already at a few hundreds of digits it is preferable to load
`irwin_v5.sage` which will try to take advantage of multiple cores on your
system.  It defaults to using `8` workers, but this can be configured by
setting `maxworkers` variable *and reloading* `irwin_v5.sage`.

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

## Additional details

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

- Here are some timings of `irwin()` from `irwin_v5` module, on a 
  Mac Mini M4 Pro with `10+4` cores:
  - version 1.5.4 of April 23, 2025, computed `2+101010` decimals of the "no-9"
    Kempner constant in `5h38mn`, with `maxworkers` left to its default `8`.
  - version 1.5.5 of April 30, 2025, computed `2+130010` decimals in `10h10mn`,
    using `maxworkers=10`.

- The next files with names of the type `k_prec_N` contain decimal expansions
  of the classic "no-9 radix-10" Kempner series `22.92067661926415...`,
  correctly rounded to `N` decimal places.  They are currently (sorry for the
  `+` in the filenames which is trying to say how many digits before and after
  decimal separator):
  * [k_prec_2+1000](k_prec_2+1000)
  * [k_prec_2+2000](k_prec_2+2000)
  * [k_prec_2+5000](k_prec_2+5000)
  * [k_prec_2+10000](k_prec_2+10000)
  * [k_prec_2+20000](k_prec_2+20000)
  * [k_prec_2+25000](k_prec_2+25000)
  * [k_prec_2+50000](k_prec_2+50000)
  * [k_prec_2+99998](k_prec_2+99998) (contributed by Yusuf Emin Akpınar, thanks!)
  
- [taille_pascal.pdf](taille_pascal.pdf) explains how many bits are needed to
  store in computer memory the Pascal triangle up (or rather down) to a
  certain row, or the memory needed to store only one such row.  Thanks to
  Nicolas Radulesco who asked for SageMath equivalents of the few Maple usages
  which were inserted in the text, they are in
  [taille_pascal_symbolic](taille_pascal_symbolic).

## The `@parallel` issue #1

- I should note here that there are some surprises with parallelization using
  `@parallel` on my very recently acquired macOS 15.4.1 Sequoia.  To
  parallelize the computation by recurrence of the coefficients which I have
  associated with Kempner and Irwin sums, the code necessarily has to call
  many times a procedure which has been `@parallel`-decorated; for example if
  we need `10000` such coefficients and use `8` workers, we will call `1250`
  times the `@parallel`-ized procedure each time with an input having `8`
  entries.

  Turns out that on macOS 15.4.1 at least (but I can't test easily on other
  systems), as I discovered in issue #1, and as is explored in test files such
  as [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage), iterating
  calls to a procedure calling a `@parallel` decorated one causes a seemingly
  linear increase of the execution time at each new call.  With
  [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage) this time drift
  seems to continue more or less forever, in more realistic situations it is
  seen to get reset from time to time (probably in connection with some "cache
  cleanup" somewhere by the OS).  The initial faster execution times will
  never be matched again though.  Only way known to author at time of
  writing is to quit the SageMath interactive session and to relaunch a fresh
  one.  Here is a typical illustration of the execution times drifting,
  where the workers only action is to "sleep"!

  ```text
  $ sage
  ┌────────────────────────────────────────────────────────────────────┐
  │ SageMath version 10.6, Release Date: 2025-03-31                    │
  │ Using Python 3.12.5. Type "help()" for help.                       │
  └────────────────────────────────────────────────────────────────────┘
  sage: load("test_parallel_sleep_v2.sage")

  The bar() function is set-up to use per default numcalls=1000
  and T=0.00100000000000000.  Its execution will time consecutively
  chunks of 100 calls to the @parallel decorated interface
  to time.sleep().

  The variable "ncpus" sets ncpus for @parallel usage. You
  can set it and then redo load("test_parallel_sleepv2.sage").

  ncpus variable has been created and assigned value 8.

  sage: %time bar()
  Using ncpus=8 parallel decorated function
  1.131s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001289s.
  1.154s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001317s.
  1.209s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001386s.
  1.251s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001439s.
  1.267s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001458s.
  1.337s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001546s.
  1.394s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001617s.
  1.427s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001659s.
  1.442s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001677s.
  1.450s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001688s.
  CPU times: user 398 ms, sys: 9.75 s, total: 10.1 s
  Wall time: 13.1 s
  sage: %time bar()
  Using ncpus=8 parallel decorated function
  1.558s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001823s.
  1.594s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001868s.
  1.630s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001913s.
  1.711s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002014s.
  1.783s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002104s.
  1.845s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002181s.
  1.855s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002194s.
  1.906s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002258s.
  1.949s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002311s.
  2.003s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002379s.
  CPU times: user 411 ms, sys: 14.3 s, total: 14.7 s
  Wall time: 17.8 s
  sage: %time bar()
  Using ncpus=8 parallel decorated function
  2.082s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002478s.
  2.049s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002436s.
  2.111s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002513s.
  2.153s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002567s.
  2.187s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002608s.
  2.429s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002911s.
  2.519s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003024s.
  2.710s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003262s.
  2.854s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003442s.
  2.980s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003600s.
  CPU times: user 424 ms, sys: 20 s, total: 20.4 s
  Wall time: 24.1 s
  sage:
  ```
  
  It looks as if this is tied with heavy IO, as the macOS Activity Monitor
  shows disk writes at about 1Mo/s (starting at more than 2Mo/s) which I guess
  is related to pickling the full Sage state...

  I have tested using the `multiprocessing` library via its `Pool.map()` and
  the starting methods `'fork'`, `'forkserver'` or `'spawn'` as well as the
  `concurrent.futures.ProcessPoolExecutor()` and I have never been able to
  reproduce any systematic drift in execution times, contrarily to what
  happens with `@parallel`.
  
  Sadly though, all of these methods, when tried out with the computation of
  the coefficients at the heart of my research, cause higher overhead. My
  algorithm in `irwin_v5.sage` tests regularly if it is worthwile to go
  parallel for the recurrences computing the coefficients defined in my
  research, and for these new protocols the decision to toggle on the parallel
  computations is never taken.  Hence we are back for that part to the same
  execution times as with the `irwin_v3.sage` more or less, and gains are
  obtained only for the later part which computes the "`beta(m+1)`'s".

  If only there was a way after each call to a `@parallel`-decorated procedure
  to tell the OS to clean up its stuff so that we go back to similar state as on
  first execution!  Probably there is a way but for this author it is quite a
  challenge to dig into these things.  It can't be a widespread issue with
  `@parallel` and must be limited probably to recent macOSes and perhaps only
  ARM architectures, else people would have reported it. (But one is surprised
  sometimes for such things; the author of these lines has more than once
  found bugs in decades old software.)

  I also tested in `C` using the `OpenMP` library if any time drift could be
  put into evidence, and none showed up.  Thanks to Yusuf Emin Akpınar for
  providing `C` implementation of the algorithm of my research papers.
  
  I wish I could test on more systems, only thing I can say is that nothing of
  this sort is seen on old macOSes such as High Sierra.  I am quite interested
  into reports on how
  [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage) behaves on
  various systems.  Please comment at #1.

  In practice, what this means is that on afflicted systems, one should quit
  `sage` and restart a session before launching a run for which one needs the
  best efficiency.  But it also means that one-shot runs could probably be
  faster because the time drift does affect the successive calls internally
  needed to build up the recurrent coefficients (cf. for example
  [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage)).

  These files are related to testing the `@parallel` issue #1:
  * [test_openmp_sleep.c](test_openmp_sleep.c)
  * [test_parallel_sleep.sage](test_parallel_sleep.sage)
  * [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage)
  * [test_executor_sleep.sage](test_executor_sleep.sage)
  * [test_pool_fork_sleep.py](test_pool_fork_sleep.py)
  * [test_pool_spawn_sleep.py](test_pool_spawn_sleep.py)
  * [test_pool_fork_RF.sage](test_pool_fork_RF.sage)
  * [test_pool_fork_RF.py](test_pool_fork_RF.py)
  * [test_pool_spawn_RF.py](test_pool_spawn_RF.py)
  * [test_pool_forkserver_RF.py](test_pool_forkserver_RF.py)

## Bibliographical references

This repository is devoted to the formulas discovered in:

- [Moments in the exact summation of the curious series of Kempner type](https://arxiv.org/abs/2402.08525)
- [Measures for the summation of Irwin series](https://arxiv.org/abs/2402.09083)

Earlier numerical works computing otherwise Kempner and Irwin sums include:

- Robert Baillie: Sums of reciprocals of integers missing a given digit. Amer. Math.
Monthly 86(5), 372–374 (1979) [DOI](https://doi.org/10.2307/2321096Sums)
- Robert Baillie: Summing the curious series of Kempner and Irwin (2008).
[arXiv:0806.4410](https://arxiv.org/abs/0806.4410)
- Thomas Schmelzer and Robert Baillie: Summing a curious, slowly convergent series. Amer.
Math. Monthly 115(6), 525–540 (2008) [DOI](https://doi.org/10.1080/00029890.2008.11920559)

My two papers quoted above have not yet been published but some of my further
research already has appeared:

- Summing the "exactly one 42" and similar subsums of the harmonic series,  Advances in Applied Mathematics Volume 162, January 2025, 102791. [DOI](https://doi.org/10.1016/j.aam.2024.102791)
- Digamma function and general Fischer series in the theory of Kempner sums, Expositiones Mathematicae, Volume 42, Issue 6, December 2024, 125604. [DOI](https://doi.org/10.1016/j.exmath.2024.125604)
- Measures associated with certain ellipsephic harmonic series and the Allouche-Hu-Morin limit theorem, Acta Mathematica Hungarica (2025) [DOI](https://doi.org/10.1007/s10474-025-01525-3)

You will find at
[arXiv:burnol](https://arxiv.org/search/?searchtype=author&query=Burnol%2C+J)
some further manuscripts related to this topic, as well as earlier papers
doing fancy mathematics on some other topics.  For fancier earlier mathematics still,
the links are not yet public.


TODO: perhaps add some additional bibliographical references and do the
fancier earlier mathematics at long last.

## License

The files in this repository are distributed under the
CC-BY-SA 4.0 License.  See [LICENSE](LiCENSE).
