@parallel(ncpus=8)
def foo(x):
    return x**2

def test(R, check=False):
    results = list(foo((a,) for a in range(R)))
    squares = [result[1] for result in sorted(results)]
    if check:
        print(squares)
    return None

