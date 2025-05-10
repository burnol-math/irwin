# -*- mode: python -*-
# load("test_parallel_spawn_sleep.sage")
# Cost of spawning is prohibitive so numcalls_default is
# set to only 1, and each call is timed...
import time

T_default = 0.001
numcalls_default = 10
ncpus_default = 8

def bar(numcalls=numcalls_default,
        T=T_default,
        ncpus=ncpus_default):

    @parallel(ncpus=ncpus, p_iter='multiprocessing')
    def foo(a, fT, time):
        time.sleep(fT)
        return None

    fT = float(T)
    print(f"Using ncpus={ncpus} parallel decorated function")
    for _ in range(numcalls):
        print(f"Now spawning @parallelized sleep({fT:.3f}) on {ncpus} workers.",
              end =" ",
              flush = True)
        start = time.perf_counter()
        results = list(foo((a, fT, time) for a in range(ncpus)))
        delta = time.perf_counter()-start
        print(f"{delta:.3f}s.\n Cost/ncpus is {(delta - T)/ncpus:.6f}s.")
        
    return None

if __name__ == "__main__":
    print(f"""
numcalls={numcalls_default}
T={T_default}
ncpus={ncpus_default}
You can set these variables as parameters to the bar() function.
The latter calls numcalls times a function foo() which sleeps
T seconds, which is @parallel decorated with p_iter set
to 'multiprocessing' and ncpus workers.  This is an inner function
to bar() so ncpus can become an argument of the latter.
So foo is applied to a list of ncpus entries, and this is done
numcalls number of times.  The default for numcalls is small
due to the overhead of spawning.  Start with bar() or %time bar()."""
          )
