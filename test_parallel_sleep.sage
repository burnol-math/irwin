# -*- mode: python -*-
# https://github.com/sagemath/sage/issues/39960#issuecomment-2839226392
ncpus = 8

import time

@parallel(ncpus=ncpus)
def foo(a, T):
    time.sleep(float(T))
    return None

def bar(numcalls=100, T=0.005):
    for i in range(numcalls):
        results = list(foo((a, T,) for a in range(ncpus)))
    return None
