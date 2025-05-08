# -*- mode: python -*-
# https://github.com/sagemath/sage/issues/39960#issuecomment-2839226392

import time

try:
    _ = ncpus
    ncpusinfostring = ("ncpus variable already existed "
                       f"with value {ncpus}.")
except NameError:
    ncpus = 8
    ncpusinfostring = ("ncpus variable has been created "
                       "and assigned value 8.")


@parallel(ncpus=ncpus)
def foo(a, fT):
    time.sleep(fT)
    return None

T_default = 0.001
numcalls_default = 1000

def bar(numcalls=numcalls_default, T=T_default):
    fT = float(T)
    print(f"Using ncpus={ncpus} parallel decorated function")
    q, r = divmod(numcalls, 100)
    for _ in range(q):
        start = time.perf_counter()
        for _ in range(100):
            results = list(foo((a, fT) for a in range(ncpus)))
        delta = time.perf_counter()-start
        print(f"{delta:.3f}s per 100 x {ncpus} @para sleep({fT:.3f}). "
              f"Cost/ncpus/calls is {(delta/100 - T)/ncpus:.6f}s.")
    if r > 0:
        start = time.perf_counter()
        for _ in range(r):
            results = list(foo((a, fT) for a in range(ncpus)))
        delta = time.perf_counter()-start
        print(f"{delta:.3f}s per {r} x {ncpus} @para sleep({fT:.3f}). "
              f"Cost/ncpus/calls is {(delta/r - T)/ncpus:.6f}s.")
        
    return None

if __name__ == "__main__":
    print(f"""
The bar() function is set-up to use per default numcalls={numcalls_default}
and T={T_default}.  Its execution will time consecutively
chunks of 100 calls to the @parallel decorated interface
to time.sleep().

The variable "ncpus" sets ncpus for @parallel usage. You
can set it and then redo load("test_parallel_sleepv2.sage").

{ncpusinfostring}
"""
          )
