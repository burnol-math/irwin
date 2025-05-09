# -*- mode: python -*-
# test_pool_fork_RF.py

# sage -ipython  # un iPython qui sait trouver les librairies
#                # de SageMath-10-6.app sur macOS
# Puis par exemple:
# from test_pool_fork_RF import bar as fork_bar

from sage.rings.real_mpfr import RealField
import time
import multiprocessing as mp

def worker(args):
    a, data, precision = args
    R = RealField(precision)
    x = R(data)
    result = x.sin()
    return result

def bar(numcalls=100,
        ncpus=8,
        data_list=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
        precision=100,
        check=False):
    mp.set_start_method('fork', force=True)
    print("Setting multiprocessing start method to 'fork'")
    with mp.Pool(processes=ncpus) as pool:
        for i in range(numcalls):
            inputs = [(a, data_list[a], precision) for a in range(ncpus)]
            results = pool.map(worker, inputs)
    if check:
        print(results)
    return None

if __name__ == "__main__":
    mp.set_start_method('fork', force=True)
    print("Setting multiprocessing start method to 'fork'")

