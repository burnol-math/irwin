# The `@parallel` issue #1

There are some surprises with parallelization using
`@parallel` on my very recently acquired macOS 15.4.1 Sequoia.  To
parallelize the computation by recurrence of the coefficients which I have
associated with Kempner and Irwin sums, the code necessarily has to call
many times a procedure which has been `@parallel`-decorated; for example if
we need `10000` such coefficients and use `8` workers, we will call `1250`
times the `@parallel`-ized procedure each time with an input having `8`
entries.

Turns out that on macOS 15.4.1 at least (but I can't test easily on other
systems), as I discovered in issue #1, and as is explored in test files such
as [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage), iterating
calls to a procedure calling a `@parallel` decorated one causes a seemingly
linear increase of the execution time at each new call.  With
[test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage) this time drift
seems to continue more or less forever, in more realistic situations it is
seen to get reset from time to time (probably in connection with some "cache
cleanup" somewhere by the OS).  The initial faster execution times will
never be matched again though.  Only way known to author at time of
writing is to quit the SageMath interactive session and to relaunch a fresh
one.  Here is a typical illustration of the execution times drifting,
where the workers only action is to "sleep"!

```text
$ sage
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 10.6, Release Date: 2025-03-31                    │
│ Using Python 3.12.5. Type "help()" for help.                       │
└────────────────────────────────────────────────────────────────────┘
sage: load("test_parallel_sleep_v2.sage")

The bar() function is set-up to use per default numcalls=1000
and T=0.00100000000000000.  Its execution will time consecutively
chunks of 100 calls to the @parallel decorated interface
to time.sleep().

The variable "ncpus" sets ncpus for @parallel usage. You
can set it and then redo load("test_parallel_sleepv2.sage").

ncpus variable has been created and assigned value 8.

sage: %time bar()
Using ncpus=8 parallel decorated function
1.131s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001289s.
1.154s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001317s.
1.209s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001386s.
1.251s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001439s.
1.267s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001458s.
1.337s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001546s.
1.394s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001617s.
1.427s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001659s.
1.442s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001677s.
1.450s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001688s.
CPU times: user 398 ms, sys: 9.75 s, total: 10.1 s
Wall time: 13.1 s
sage: %time bar()
Using ncpus=8 parallel decorated function
1.558s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001823s.
1.594s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001868s.
1.630s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001913s.
1.711s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002014s.
1.783s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002104s.
1.845s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002181s.
1.855s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002194s.
1.906s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002258s.
1.949s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002311s.
2.003s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002379s.
CPU times: user 411 ms, sys: 14.3 s, total: 14.7 s
Wall time: 17.8 s
sage: %time bar()
Using ncpus=8 parallel decorated function
2.082s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002478s.
2.049s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002436s.
2.111s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002513s.
2.153s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002567s.
2.187s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002608s.
2.429s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002911s.
2.519s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003024s.
2.710s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003262s.
2.854s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003442s.
2.980s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003600s.
CPU times: user 424 ms, sys: 20 s, total: 20.4 s
Wall time: 24.1 s
sage:
```

One observes during execution heavy disk write activity, as reported by the
macOS Activity Monitor at about 1Mo/s throughout, starting at more than
2Mo/s.  But this may be not really relevant as the pure Python reproducer
(which is commented below) does not trigger any disk activity but does
show the timings drift.

I have tested using the `multiprocessing` library via its `Pool.map()` and
the starting methods `'fork'`, `'forkserver'` or `'spawn'` as well as the
`concurrent.futures.ProcessPoolExecutor()` and I have never been able to
reproduce any systematic drift in execution times, contrarily to what
happens with `@parallel`.

Sadly though, all of these methods, when tried out with the computation of
the coefficients at the heart of my research, cause higher overhead. My
algorithm in `irwin_v5.sage` tests regularly if it is worthwile to go
parallel for the recurrences computing the coefficients defined in my
research, and for these new protocols the decision to toggle on the parallel
computations is never taken.  Hence we are back for that part to the same
execution times as with the `irwin_v3.sage` more or less, and gains are
obtained only for the later part which computes the "`beta(m+1)`'s".

If only there was a way after each call to a `@parallel`-decorated procedure
to tell the OS to clean up its stuff so that we go back to similar state as on
first execution!  Probably there is a way but for this author it is quite a
challenge to dig into these things.  It can't be a widespread issue with
`@parallel` and must be limited probably to recent macOSes and perhaps only
ARM architectures, else people would have reported it. (But one is surprised
sometimes for such things; the author of these lines has more than once
found bugs in decades old software.)

I also tested in `C` using the `OpenMP` library if any time drift could be
put into evidence, and none showed up.  Thanks to Yusuf Emin Akpınar for
providing `C` implementation of the algorithm of my research papers.

I wish I could test on more systems, only thing I can say is that nothing of
this sort is seen on old macOSes such as High Sierra.  I am quite interested
into reports on how
[test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage) behaves on
various systems.  Please comment at #1.

In practice, what this means is that on afflicted systems, one should quit
`sage` and restart a session before launching a run for which one needs the
best efficiency.  But it also means that one-shot runs could probably be
faster because the time drift does affect the successive calls internally
needed to build up the recurrent coefficients (cf. for example
[test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage)).

The file [test_pyparallel_sleep.py](test_pyparallel_sleep.py) is a
reproducer in pure Python of this problem.  Thanks to Edgar Costa for
providing the `@parallel` emulation in pure Python which it uses.
Here is with macOS 15.5 and Python 3.13:

```text
$ ipython
Python 3.13.3 (v3.13.3:6280bb54784, Apr  8 2025, 10:47:54) [Clang 15.0.0 (clang-1500.3.9.4)]
Type 'copyright', 'credits' or 'license' for more information
IPython 9.2.0 -- An enhanced Interactive Python. Type '?' for help.
Tip: You can change the editing mode of IPython to behave more like vi, or emacs.

In [1]: from test_pyparallel_sleep import bar

In [2]: for _ in range(5):
   ...:     %time bar()
   ...: 
Using ncpus=8 parallel decorated function
0.675s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.000718s.
0.738s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.000798s.
0.799s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.000873s.
0.891s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.000988s.
0.945s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001057s.
1.018s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001148s.
1.082s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001227s.
1.152s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001314s.
1.233s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001416s.
1.312s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001516s.
CPU times: user 226 ms, sys: 6.92 s, total: 7.14 s
Wall time: 9.85 s
Using ncpus=8 parallel decorated function
1.405s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001632s.
1.480s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001725s.
1.562s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001828s.
1.628s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001910s.
1.698s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.001997s.
1.775s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002093s.
1.860s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002200s.
1.961s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002326s.
2.057s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002446s.
2.140s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002549s.
CPU times: user 222 ms, sys: 14 s, total: 14.2 s
Wall time: 17.6 s
Using ncpus=8 parallel decorated function
2.227s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002658s.
2.316s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002770s.
2.403s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002879s.
2.478s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.002973s.
2.557s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003072s.
2.654s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003192s.
2.754s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003318s.
2.852s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003440s.
2.920s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003525s.
3.027s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003659s.
CPU times: user 238 ms, sys: 21.4 s, total: 21.7 s
Wall time: 26.2 s
Using ncpus=8 parallel decorated function
3.133s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003791s.
3.225s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.003906s.
3.312s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004015s.
3.391s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004113s.
3.496s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004245s.
3.565s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004332s.
3.744s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004555s.
4.029s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.004912s.
4.563s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005578s.
4.312s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005265s.
CPU times: user 281 ms, sys: 30.7 s, total: 31 s
Wall time: 36.8 s
Using ncpus=8 parallel decorated function
4.473s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005466s.
4.355s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005319s.
4.410s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005388s.
4.656s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005695s.
4.629s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005662s.
4.530s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005537s.
4.705s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.005756s.
4.900s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.006000s.
4.932s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.006040s.
5.034s per 100 x 8 @para sleep(0.001). Cost/ncpus/calls is 0.006167s.
CPU times: user 310 ms, sys: 39.5 s, total: 39.8 s
Wall time: 46.6 s
```

These further files are also related to testing the `@parallel` issue #1:
* [test_openmp_sleep.c](test_openmp_sleep.c)
* [test_parallel_sleep.sage](test_parallel_sleep.sage)
* [test_parallel_sleep_v2.sage](test_parallel_sleep_v2.sage)
* [test_executor_sleep.sage](test_executor_sleep.sage)
* [test_pool_fork_sleep.py](test_pool_fork_sleep.py)
* [test_pool_spawn_sleep.py](test_pool_spawn_sleep.py)
* [test_pool_fork_RF.sage](test_pool_fork_RF.sage)
* [test_pool_fork_RF.py](test_pool_fork_RF.py)
* [test_pool_spawn_RF.py](test_pool_spawn_RF.py)
* [test_pool_forkserver_RF.py](test_pool_forkserver_RF.py)
