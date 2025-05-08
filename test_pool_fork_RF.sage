# -*- mode:python -*-
# test_parallel_mp_pool_fork.sage
import time
import multiprocessing as mp

mp.set_start_method('fork', force=True)

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
    with mp.Pool(processes=ncpus) as pool:
        for i in range(numcalls):
            inputs = [(a, data_list[a], precision) for a in range(ncpus)]
            results = pool.map(worker, inputs)
    if check:
        print(results)
    return None
