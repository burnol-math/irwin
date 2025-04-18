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

# Trying to be somewhat close to real use case in irwin_v4
@parallel(ncpus=ncpus)
def machin(x):
    return [1/RealField(z[1] - 3*z[0])(999 ** (z[0]+1)) for z in x]

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
def test5(R, check=False, fn=machin):
    if fn == machin:
        Pmax = 3 * R + 20
        inputs = list(zip(range(R), [Pmax] * R ))
    else:
        inputs = list(range(R))
    blocks = [inputs[i::ncpus] for i in range(ncpus)]
    results_1 = [result[1] for result in sorted(list(fn(blocks)))]
    delta = len(results_1[0])-len(results_1[-1])
    if delta > 0:
        results_1[-1].extend([None] * delta)
        powers = [x for xs in zip(*results_1) for x in xs][:-delta]
    else:
        powers = [x for xs in zip(*results_1) for x in xs]
    if check:
        print(powers)
    # do some (very) minimal check always
    # powers[0] occupies roughly of the order of R decimal
    # digits if fn==machin
    # Attention que le type double donné par float()
    # a en fait un range assez faible d'exposants scientifiques
    # possibles... 1/float(10**324) donne 0.0
    return float(powers[10]) if R>10 else powers[0], powers[-1]

    
