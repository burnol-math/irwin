#!/usr/bin/env python3
# test_pyparallel_sleep.py
# ipython -i test_pyparallel_sleep.py
# or
# python -i test_pyparallel_sleep.py
# or
# from test_pyparallel_sleep import bar
"""
NOTA BENE: THE CODE FOR THE parallel DECORATOR IS ENTIRELY TAKEN
FROM
# 29 avril 2025
# https://gist.githubusercontent.com/edgarcosta/77d09c3727da202153cbcf93cc828177/raw/1597883c15b49f732771dc511b82144fc3fc3e7d/parallel_purepython.py
THANKS!

HERE IS ITS ORIGINAL DOCSTRING:
Minimal fork-based parallel decorator.
• Unix-only (needs os.fork).
• Each job is a tuple of positional arguments.
• When you call the decorated function with an *iterable* of
  those tuples it returns an *iterator* of results (in the order
  each child finishes).

The rest of the code is a mix of these files from this repo:

    test_parallel_sleep_v2.sage
    test_pool_fork_sleep.py

The main change is that ncpus was made an argument to bar()
and the decorated foo() is made an inner function.  I hope
this means decoration will use the ncpus as in bar arglist.
"""
import os, pickle, signal, time
from types import GeneratorType
from typing import Iterable, Iterator, Tuple, Any

# --------------------------------------------------------------------------- #
# decorator                                                                   #
# --------------------------------------------------------------------------- #
def parallel(ncpus: int | None = None):
    """@parallel(ncpus=4) decorator factory."""

    ncpus = ncpus or os.cpu_count()

    def _decorator(func):
        def _runner(arg_iter: Iterable[Tuple[Any, ...]]) -> Iterator[Any]:
            arg_iter = iter(arg_iter)          # so we can call next()

            # pid  -> read-fd  (one pipe per live child)
            active: dict[int, int] = {}

            while True:
                # ---------- spawn up to ncpus workers --------------------- #
                while len(active) < ncpus:
                    try:
                        args = next(arg_iter)
                    except StopIteration:
                        break                   # no more work
                    r_fd, w_fd = os.pipe()
                    pid = os.fork()

                    if pid == 0:                # ── CHILD ─────────────────
                        os.close(r_fd)
                        try:
                            rv = func(*args)
                            payload = pickle.dumps((True, rv),
                                                    protocol=pickle.HIGHEST_PROTOCOL)
                        except Exception as e:   # send exception back
                            import traceback
                            tb = traceback.format_exc()
                            payload = pickle.dumps((False, (e, tb)),
                                                    protocol=pickle.HIGHEST_PROTOCOL)
                        os.write(w_fd, payload)
                        os.close(w_fd)
                        os._exit(0)

                    # ── PARENT ─────────────────────────────────────────────
                    os.close(w_fd)
                    active[pid] = r_fd

                if not active:                  # no running or queued jobs
                    break

                # ---------- reap one finished child ----------------------- #
                done_pid, _ = os.wait()         # waits for *any* child
                r_fd = active.pop(done_pid)

                # read whatever the child wrote
                data = b""
                while chunk := os.read(r_fd, 4096):
                    data += chunk
                os.close(r_fd)

                ok, payload = pickle.loads(data)
                if ok:
                    yield payload               # normal result
                else:                           # re-raise child exception
                    exc, tb = payload
                    raise RuntimeError(f"child process raised:\n{tb}") from exc

        return _runner
    return _decorator

# --------------------------------------------------------------------------- #
# From a mix of test_parallel_sleep_v2.sage with test_pool_fork_sleep.py
# The decorated function has been made an inner function so that it is
# easier to modify ncpus in the decoration.
# --------------------------------------------------------------------------- #

T_default = 0.001
numcalls_default = 1000
ncpus_default = 8

def bar(numcalls=numcalls_default,
        T=T_default,
        ncpus=ncpus_default,
        check=False):
    # make this an inner function to pass ncpus as a variable
    @parallel(ncpus=ncpus)  # <<< using the pure python decorator
    def foo(a, fT, check):
        time.sleep(fT)
        return None if not check else a**2

    # fT = float(T)
    fT = T
    print(f"Using ncpus={ncpus} parallel decorated function")
    q, r = divmod(numcalls, 100)
    m = 0
    for _ in range(q):
        start = time.perf_counter()
        for _ in range(100):
            results = list(foo((a, fT, check) for a in range(m, m+ncpus)))
        delta = time.perf_counter()-start
        print(f"{delta:.3f}s per 100 x {ncpus} @para sleep({fT:.3f}). "
              f"Cost/ncpus/calls is {(delta/100 - T)/ncpus:.6f}s.")
        m += ncpus
    if r > 0:
        start = time.perf_counter()
        for _ in range(r):
            results = list(foo((a, fT, check) for a in range(m, m+ncpus)))
        delta = time.perf_counter()-start
        print(f"{delta:.3f}s per {r} x {ncpus} @para sleep({fT:.3f}). "
              f"Cost/ncpus/calls is {(delta/r - T)/ncpus:.6f}s.")
        m += ncpus
    if check:
        # Only print results for last batch.
        # In brief testing the results come out already in order as in the
        # inputs, which is probably to be expected as per parallel comments.
        print(results)

    return None

if __name__ == "__main__":
    print(f"""
numcalls={numcalls_default}
T={T_default}
ncpus={ncpus_default}
You can set these variables as parameters to bar().
With iPython, you can do this:

for _ in range(3):
    %time bar()

With python, drop the %time part, anyhow bar()
furnishes its own time measurements."""
          )
