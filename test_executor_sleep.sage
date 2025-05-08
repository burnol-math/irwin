# -*- mode: python -*-
import multiprocessing
import concurrent.futures
import time

try:
    _ = maxworkers
    maxworkersinfostring = ("maxworkers variable already existed "
                       f"with value {maxworkers}.")
except NameError:
    maxworkers = 8
    maxworkersinfostring = ("maxworkers variable has been created "
                       "and assigned value 8.")

def foo(args):
    a, fT, check = args
    time.sleep(fT)
    # hopefully a**2 is not computed if check is False
    return None if not check else a**2

T_default = 0.001
numcalls_default = 1000

def bar(numcalls=numcalls_default, T=T_default, maxworkers=maxworkers, check=False):
    fT = float(T)
    print(f"Using maxworkers={maxworkers} as max_workers for ProcessPoolExecutor")
    with concurrent.futures.ProcessPoolExecutor(
            mp_context=multiprocessing.get_context("fork"),
            max_workers=maxworkers,
    ) as executor:
        m = 0
        q, r = divmod(numcalls, 100)
        for j in range(q):
            start = time.perf_counter()
            for i in range(100):
                results = list(executor.map(foo, [(a,
                                                   fT,
                                                   check) for a in range(m,m+maxworkers)]))
                m += maxworkers
            delta = time.perf_counter()-start
            print(f"{delta:.3f}s per 100 x {maxworkers} @para sleep({fT:.3f}). "
                  f"Cost/maxworkers/calls is {(delta/100 - T)/maxworkers:.6f}s.")
            if check:
                print(results)
        if r > 0:
            start = time.perf_counter()
            for i in range(r):
                results = list(executor.map(foo, [(a,
                                                   fT,
                                                   check) for a in range(m,m+maxworkers)]))
                m += maxworkers
            delta = time.perf_counter()-start
            print(f"{delta:.3f}s per {r} x {maxworkers} @para sleep({fT:.3f}). "
                  f"Cost/maxworkers/calls is {(delta/r - T)/maxworkers:.6f}s.")
            if check:
                print(results)
       
    return None

if __name__ == "__main__":
    print(f"""
The bar() function is set-up to use per default numcalls={numcalls_default}
and T={T_default}.  Its execution will time consecutively
chunks of 100 calls each mapping function foo() (which sleeps)
via concurrent.futures.ProcessPoolExecutor() on a list of
exactly {maxworkers} (dummy, for this test) arguments.

{maxworkersinfostring}
"""
          )
