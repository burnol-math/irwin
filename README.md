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
  [irwin_loader.py](irwin_loader.py) and similarly named file can be used to
  load concurrently `irwin.sage` and `irwin_v3.sage` into the same interactive
  SageMath session.  See the comments therein.  No guarantees though.

- `irwin_v4.sage` (Not Yet) adds some "parallelization" attempt.  This is
  experimental and mainly effective when computing the `beta(m+1)`'s.  *attow
  this file is not yet pushed to the repository*.

- The next files with names of the type `k_prec_N` contain decimal expansions
  of the classic "no-9 radix-10" Kempner series `22.92067661926415...`,
  correctly rounded to `N` decimal places.  They are currently:
  * [k_prec_1000](k_prec_1000)
  * [k_prec_10000](k_prec_10000)
  * [k_prec_20000](k_prec_20000)
  * [k_prec_25000](k_prec_25000)
  * [k_prec_50000](k_prec_50000)
  * [k_prec_99998](k_prec_99998) (contributed by Yusuf Emin Akpınar, thanks!)
  
- [taille_pascal.pdf](taille_pascal.pdf) explains how many bits are needed to
  store in computer memory the Pascal triangle up (or rather down) to a
  certain row, or the memory needed to store only one such row.

## License

The contributed files (Python, SageMath) are distributed under the
CC-BY-SA 4.0 License.  See [LICENSE](LiCENSE).
