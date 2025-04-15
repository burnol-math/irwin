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
banner instructions which will be printed.

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
  needing binomial coefficients up to `N=25000`) circa 470 gigabytes of
  computer memory.  In contrast, the 2025 version needs less than 113
  megabytes...

- `irwin_v4.sage` (Not Yet) adds some "parallelization" attempt.  This is
  experimental and mainly effective when computing the `beta(m+1)` sums of
  inverse powers on integers having a number of digits equal to the `level`
  parameter of the algorithm.  *attow this file is not yet pushed to the
  repository*.


## License

The contributed files (Python, SageMath) are distributed under the
CC-BY-SA 4.0 License.  See [LICENSE](LiCENSE).
