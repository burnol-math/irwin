# -*- mode: python -*-

# NOTA BENE (May 1, 2025)

# This was done in April at some point in development of
# irwin_v4.sage to test relative efficiencies of using
# @parallel.  For some forgotten rasons, the test2(), test3(),
# test4() and test5() all call the @parallel function with some
# explicit list range, rather than simply a start, end, step
# data. This probably causes some unefficiency and anyhow
# irwin_v5.sage has improved on this later on.  But this is not
# backported here.

# Another very important point to make is that there is
# definitely an upstream bug somewhere which causes usage or
# @parallel decorated functions to have very counter-intuitive
# after effects -- **after** the calls are done... See issue #1.
# This may interfere with various benchmarking which was
# attempted with this test file.

# Finally the procedure machin() does not emulate good enough
# the real life things done by irwin_v4 and later.

ncpus=8

@parallel(ncpus=ncpus)
def foo(x):
    return x**2

@parallel(ncpus=ncpus)
def bar2(x):
    return [z**2 for z in x]

@parallel(ncpus=ncpus)
def bar100(x):
    # use something which should be more costly for some input
    # not enough practitioner of Python to make an informed
    # choice here
    return [z**100 for z in x]

# Trying to be somewhat close to real use case in irwin_v4. But
# in powers of 999 was perhaps not a good idea and have been
# replaced by powers of say 101.
@parallel(ncpus=ncpus)
def machin(x):
    return [1/RealField(z[1] - 3*z[0])(101 ** (z[0]+1)) for z in x]

# BAD WAY
def test(R, check=False):
    results = list(foo((a,) for a in range(R)))
    squares = [result[1] for result in sorted(results)]
    if check:
        print(squares)
    return None

# CORRECT WAY
# https://github.com/sagemath/sage/issues/39960#issuecomment-2814217741
def test2(R, check=False):
    inputs = list(range(R))
    blocks = [inputs[i::ncpus] for i in range(ncpus)]
    results = list(bar2(blocks))
    squares = sum([result[1] for result in sorted(results)], [])
    if check:
        print(squares)
    return None

# But I need the output ordered in exactly same as input,
# hence the following test3(), test4(), test5().
#
# When only the length R matters, it seems that the best
# choice is test4, then test5, then test3, say with
# fn=bar2 or fn=bar100 and R tested up to 10,000,000
#
# But using machin for fn is closer to what matters for irwin
# and I did not measure then any significant differences.
# In particular the worry about trying to have the workload
# evenly distributed among workers as in test3() and test5()
# seems not too warranted, test4() is fine in practice
# and has simpler shorter coding.
#
def test3(R, check=False, fn=machin):
    if fn == machin:
        Pmax = 3 * R + 20
        inputs = list(zip(range(R), [Pmax] * R ))
    else:
        inputs = list(range(R))
    blocks = [inputs[i::ncpus] for i in range(ncpus)]
    results_sorted = sorted(list(fn(blocks)))
    powers = []
    q, r = divmod(R, ncpus)
    for n in range(q):
        powers.extend([result[1][n] for result in results_sorted])
    if r > 0:
        powers.extend([result[1][-1] for result in results_sorted[:r]])
    if check:
        print(powers)
        assert R == len(powers)
    # do some minimal check always
    # powers[0] occupies roughly of the order of R decimal
    # digits if fn==machin
    return float(powers[10]) if R>10 else powers[0], powers[-1]

# This way was found more efficient than test3 when raising
# integers to power 100 (fn=bar100), with R=10000000.
#
# Using machin for fn is closer to real-life and for R
# in 1000, 10000, 20000, I could not measure significant
# and obvious differences.
#
# Update May 1st: I should have tested with a procedure much
# closer to the computation of the real beta()'s because this
# testing underestimated the disadvantage of having among the 8
# parallel jobs some taking significantly longer than others.
# It is also very time-consuming and irritating to test things
# due to upstream bug/feature #1 which afflicts my system and
# makes my life lonely and miserable.  But the main problem here
# is that I did not test the real thing and for example I am
# using (this is bad) explicit lists rather than generators,
# contrarily to irwin_v5 code (which this predated anyhow).
def test4(R, check=False, fn=machin):
    q, r = divmod(R, ncpus)
    I = 0
    A = []
    for i in range(r):
        A.append(I)
        I += q + 1
    for i in range(ncpus -r):
        A.append(I)
        I += q
    assert I == R, "check your math"
    A.append(I)
    if fn == machin:
        Pmax = 3 * R + 20
        inputs = list(zip(range(R), [Pmax] * R ))
    else:
        inputs = list(range(R))
    blocks = [inputs[A[i]:A[i+1]] for i in range(ncpus)]
    results_sorted = sorted(list(fn(blocks)))
    powers = sum([result[1] for result in results_sorted], [])
    if check:
        print(powers)
        assert R == len(powers)
    # do some minimal check always
    # powers[0] occupies roughly of the order of R decimal
    # digits if fn==machin
    return float(powers[10]) if R>10 else powers[0], powers[-1]


# for fn=bar100 and R=100000 this is more efficient than test3
# but less than test4. Same for R=1,000,000 with test5 and test4
# closer. Also for R=10,000,000 test4 is the winner but it
# seems test5 is close.
#
# Using machin for fn is closer to real-life and for R
# in 1000, 10000, 20000, I could not measure significant
# and obvious differences.
#
# Sadly this test5() when I used it while working on irwin_v4
# was awfully wrong.  Fixed only much later when polishing
# irwin_v5 towards 1.5.6.
def test5(R, check=False, fn=machin):
    if fn == machin:
        Pmax = 3 * R + 20
        inputs = list(zip(range(R), [Pmax] * R ))
    else:
        inputs = list(range(R))
    blocks = [inputs[i::ncpus] for i in range(ncpus)]
    # print(blocks)
    results_1 = [result[1] for result in sorted(list(fn(blocks)))]
    # print(results_1)
    r = R % ncpus
    if r > 0:
        for j in range(1,ncpus - r + 1):
            results_1[-j].append(None)
        # print(results_1)
        # print(list(zip(*results_1)))
        powers = [x for xs in zip(*results_1) for x in xs][:-(ncpus -r)]
    else:
        powers = [x for xs in zip(*results_1) for x in xs]
    if check:
        print(powers)
        assert R == len(powers)
    # do some (very) minimal check always
    # powers[0] occupies roughly of the order of R decimal
    # digits if fn==machin
    # Attention que le type double donné par float()
    # a en fait un range assez faible d'exposants scientifiques
    # possibles... 1/float(10**324) donne 0.0
    return float(powers[10]) if R>10 else powers[0], powers[-1]
