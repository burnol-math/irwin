# test_parallel_hope.sage
import time
import multiprocessing as mp

mp.set_start_method('fork', force=True)

def worker(args):
    a, data, precision = args
    R = RealField(precision)
    x = R(data)
    result = x.sin() + x.cos()  # Example computation
    return result

ncpus = 8

def bar(numcalls=100, data_list=[0.1]*8, precision=128, ncpus=8, check=False):
    with mp.Pool(processes=ncpus) as pool:
        for i in range(numcalls):
            inputs = [(a, data_list[a], precision) for a in range(ncpus)]
            results = pool.map(worker, inputs)
    if check:
        print(results)
    return None