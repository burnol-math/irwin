ncpus=8

@parallel(ncpus=ncpus)
def foo(x):
    return x**2

@parallel(ncpus=ncpus)
def bar(x):
     return [z**2 for z in x]

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
    results = list(bar(blocks))
    squares = sum([result[1] for result in sorted(results)], [])
    if check:
        print(squares)
    return None

# DETAILS FOR MY USE CASE
def test3(R, check=False):
    inputs = list(range(R))
    blocks = [inputs[i::ncpus] for i in range(ncpus)]
    results_sorted = sorted(list(bar(blocks)))
    squares = []
    q, r = divmod(R, ncpus)
    for n in range(q):
        squares.extend([result[1][n] for result in results_sorted])
    if r > 0:
        squares.extend([result[1][-1] for result in results_sorted[:r]])
    if check:
        print(squares)
    return None

def test4(R, check=False):
    inputs = list(range(R))
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
    blocks = [inputs[A[i]:A[i+1]] for i in range(ncpus)]
    results_sorted = sorted(list(bar(blocks)))
    squares = sum([result[1] for result in results_sorted], [])
    if check:
        print(squares)
    return None

