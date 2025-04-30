# irwin


## Description

This repository is to document implementations of the formulas from the
February 2024 article
[Measures for the summation of Irwin series](https://arxiv.org/abs/2402.09083)
by me (aka Jean-François Burnol).

## Files

The list may not be always up-to-date.  So far all contributed files with
extension `.sage` are to be used interactively in a SageMath session.  Start a
SageMath session and at the `sage` prompt, enter `load("<FOO>.sage")`.  Follow
banner instructions which will be printed.  For example:

```sage
sage: load("irwin_v3.sage")
[some info printed]
sage: irwin(10,9,0,52)
22.92067661926415034816365709437593191494476243699848
```

Thanks to Arnaud Bodin and Yusuf Emin Akpınar for sharing their interest into
the Irwin paper and its predecessor
[Moments in the exact summation of the curious series of Kempner type](https://arxiv.org/abs/2402.08525),
and in particular for reporting progress in computing tens of thousands of
digits (Y. E. A. has obtained 100000=2+99998 digits of the classical "no-9"
radix 10 Kempner series).  This motivated me to revisit the 2024 code with the
perspective to make it easier and less demanding on the hardware to obtain
10000+ digits even with some high-level front-end such as SageMath.


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
  cores on the user system.  We tested `maxworkers=2` on an old hardware with
  only 2 cores, and it gave comparable result to keeping `maxworkers=8`.  The
  author is attow away from more efficient hardware and still has to test `v5`
  on such.  Keep in mind any change to `maxworkers` must be followed by a
  re-load of `irwin_v5.sage` in the interactive session.

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
