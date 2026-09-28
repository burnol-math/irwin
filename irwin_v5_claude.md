## September 28, 2026 -- Claude Opus 5.5 answers

`irwin_v5_mpmath_spawn.py` is ready. It needs `pip install mpmath gmpy2`; mpmath then uses GMP for its big integers, which matters a lot for speed. I tested it under Python 3.13.13 and 3.12, on a Linux sandbox with a single core. So I could check that results are correct, but I could not measure any parallel speed-up.

**Usage**
```python
import irwin_v5_mpmath_spawn
irwin_v5_mpmath_spawn.maxworkers = 10      # optional, taken into account at next call, no reload
from irwin_v5_mpmath_spawn import irwin, irwinpos
irwin(10, 9, 0, 1002)
```
Inside a script, the calls must sit under `if __name__ == "__main__":`. This is because the workers are started with `spawn`, the macOS default. Interactive python and ipython need nothing special.

> With the `spawn` start method, each worker is a brand-new Python interpreter, not a copy of the running one. Before the worker can do anything, `multiprocessing` makes it import your main script, under the module name `__mp_main__` instead of `__main__`. It does this so that anything defined in the script can be reconstructed in the worker.
> 
> So any code at the top level of your script runs again in every worker. Suppose the script is:
> 
> ```python
> from irwin_v5_mpmath_spawn import irwin
> print(irwin(10, 9, 0, 1002))
> ```
> 
> Here is what happens:
> 
> 1. The main process calls `irwin()`, which starts the workers.
> 2. Each worker re-imports the script, so it also executes `irwin(10, 9, 0, 1002)`.
> 3. That call in the worker tries to start its own workers.
> 4. `multiprocessing` detects this and raises a `RuntimeError`, saying a new process was started before the current one finished its bootstrapping phase. Without that safeguard, every worker would spawn workers of its own, without end.
> 
> The guard fixes this:
> 
> ```python
> from irwin_v5_mpmath_spawn import irwin
> 
> if __name__ == "__main__":
>     print(irwin(10, 9, 0, 1002))
> ```
> 
> In the main process `__name__` is `"__main__"`, so the block runs. In each worker `__name__` is `"__mp_main__"`, so the block is skipped. The worker only gets the import of `irwin_v5_mpmath_spawn`, and then waits for its tasks.
> 
> In an interactive `python` or `ipython` session there is no script file to re-import, so the workers only import `irwin_v5_mpmath_spawn`. That module's top level only defines functions and variables, and its banner is itself under a `__main__` guard. So nothing gets re-executed.
> 
> With `fork` this problem does not exist, because the child is a copy of the parent's memory and re-imports nothing. I did not use `fork` for two reasons. On macOS, Python considers it unsafe and no longer uses it as the default. Also, your pure-Python fork reproducer (`test_pyparallel_sleep.py`) is exactly the one showing the timing drift.
> 
> One side effect is worth knowing. Any expensive top-level code in your script that sits outside the guard, such as reading a large file, runs once more in each worker.


**What was converted, and how**
- **Precision.** Each `RealField(p)` became a precision `p`. Numbers are raw `mpmath.libmp` values, and every operation rounds to nearest at the precision the Sage code would have used, in the same order.
- **Integer quantities.** The `floor`/`ceil` of `nbbits_final`, `_Mmax`, `NbOfPrec` and the `IndexToR` boundaries are computed exactly with integer arithmetic.
- **Return value.** The functions return an `IrwinNumber`, which is an mpmath `mpf` rounded to `nbbits_final` bits. It prints with `nbdigits` digits, keeping trailing zeros as Sage does.
- **Unchanged parts.**
  - The signatures, the adaptive serial-versus-parallel test `(M-400)%500 < maxworkers`, and `persistentpara`.
  - The `showtimes`/`verbose` outputs and the French messages.
  - The recurrence of the $u_{j;m}$ itself.
  - The main process still keeps its Pascal rows exactly as in v5.

**Parallelization**
- **Persistent workers.** `@parallel` is replaced by `maxworkers` processes that stay alive between batches, so there is no fork at each call.
- **Recurrence.** At each parallel batch, worker $a$ receives only the rows of coefficients it does not have yet (pickled once), not the whole table. It keeps its own Pascal row, gammas and powers of $d$, so no Pascal rows are sent. This is the difference from your `irwin_v6dev_mp.sage`, which passes the whole coefficient table and a Pascal row to `pool.map` at every call.
- **Worker memory.** Each worker stores its copies rounded to the precision needed for the current $m$, and rounds them down further when the precision drops.
- **Betas.** They are split modulo `maxworkers` as in v5.

**Tests performed**
- **Reference values.** Output agrees with the values in the README and docs: 52-digit results for $k=0$, $k=1$ with all $d$, and `all=True` with $k=4$ gives the same five lines as the Sage docstring. It also matches `k_prec_2+1000` and `k_prec_2+2000` exactly, with `irwin` at levels 2 and 3 and with `irwinpos`.
- **Parallel path.** I made a test copy that forces the parallel branch at every batch. It gives bit-identical results to the serial computation for several bases, digits, values of $k$ and levels, and `irwin` agrees with `irwinpos` in all those cases.

**Worker memory for very long runs**
Because the workers are spawned, they cannot share the coefficient table copy-on-write as Sage's forks did. Each worker holds its own copy, trimmed to the current precision. For runs of 100 000 digits or more, that extra memory grows with `maxworkers`.

