# -*- mode: python -*-
# test_pool_spawn_sleep.py
# ipython -i test_pool_spawn_sleep.py
# ou pour comparer
# from test_pool_spawn_sleep import bar as spawn_bar

import time
import multiprocessing as mp

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
    return None if not check else a**2

T_default = 0.001
numcalls_default = 1000

def bar(numcalls=numcalls_default, T=T_default, check=False):
    mp.set_start_method('spawn', force=True)
    print("Setting multiprocessing start method to 'spawn'")
    # fT = float(T)
    fT = T
    print(f"Using maxworkers={maxworkers} parallel decorated function")
    q, r = divmod(numcalls, 100)
    m = 0
    with mp.Pool(processes=maxworkers) as pool:
        for _ in range(q):
            start = time.perf_counter()
            for _ in range(100):
                results = pool.map(foo, [(a, fT, check) for a in range(m, m+maxworkers)])
            delta = time.perf_counter()-start
            print(f"{delta:.3f}s per 100 x {maxworkers} sleep({fT:.3f}). "
                  f"Cost/maxworkers/calls is {(delta/100 - T)/maxworkers:.6f}s.")
            m += maxworkers
        if r > 0:
            start = time.perf_counter()
            for _ in range(r):
                results = pool.map(foo, [(a, fT, check) for a in range(m, m+maxworkers)])
            delta = time.perf_counter()-start
            print(f"{delta:.3f}s per {r} x {maxworkers} sleep({fT:.3f}). "
                  f"Cost/maxworkers/calls is {(delta/r - T)/maxworkers:.6f}s.")
            m += maxworkers
    if check:
        print(results)
        
    return None

if __name__ == "__main__":
    print(f"""
The bar() function is set-up to use per default numcalls={numcalls_default}
and T={T_default}.  Its execution will time consecutively
chunks of 100 calls to function foo() doing sleep(T) via pool.map
and starting method spawn. 

You can set maxworkers.

{maxworkersinfostring}"""
          )
    mp.set_start_method('spawn', force=True)
    print("Setting multiprocessing start method to 'spawn'")
