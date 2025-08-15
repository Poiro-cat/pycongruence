def extend_gcd(m:int, n:int):
    r0, x0, y0 = n, 0, 1
    r1, x1, y1 = m%n, 1, -(m//n)
    while r1 > 0:
        r0, r1, k = r1, r0%r1, r0//r1
        x0, x1 = x1, x0 - k*x1
        y0, y1 = y1, y0 - k*y1
    return r0, x0, y0

def power_mod(a:int, n:int, N:int):
    assert n>=0, 'Power degree n should be NOT NEGATIVE'
    assert N>0, 'Modulus N should be POSITIVE'
    y = 1
    while True:
        if n % 2 > 0: y = y*a % N
        n //= 2
        if n == 0: break
        a = a*a % N
    return y

def power_mod_quadratic_field(a:int, b:int, D:int, n:int, N:int):
    assert n>=0, 'Power degree n should be NOT NEGATIVE'
    assert N>0, 'Modulus integer N should be POSITIVE'
    u,v = 1,0
    while True:
        if n % 2 > 0: u,v = (u*a+v*b*D)%N, (u*b+v*a)%N
        n //= 2
        if n == 0: break
        a,b = (a*a+b*b*D)%N, (2*a*b)%N
    return u,v

def Legendre(a:int, p:int):
    # Please ensure input p is odd prime.
    if a % p == 0: return 0
    gcd = extend_gcd(a,p)[0]
    assert gcd == 1, f'Input p is not prime (factor = {gcd}).'
    y = power_mod(a, (p-1)//2, p)
    assert y in {1, p-1}, 'Modulus p is not prime (Euler test failed).'
    return 1 if y==1 else -1

def trial_div_factorize(N:int):
    assert N > 0, 'Factorized integer should be POSITIVE.'
    if N == 1: return set()
    factorization = {}
    y,p = N,2
    while y > 1 and p*p <= y:
        m = 0
        while y % p == 0:
            y //= p
            m += 1
        if m > 0: factorization[p] = m
        p += 1 + p%2
    if y > 1: factorization[y] = 1
    return factorization

def Sunzi(a:list[int], q:list[int]):
    assert len(a) == len(q)
    x, Q = a[0], q[0]
    for ai, qi in list(zip(a[1:],q[1:])):
        r,u,v = extend_gcd(Q,qi)
        assert r == 1, 'q not coprime'
        x = (x*qi*v + ai*Q*u) % (Q*qi)
        Q *= qi
    return x,Q

def powrt_int(n:int, k:int):
    if n == 1: return 1
    a,b = 1,n
    while b - a > 1:
        c = (a+b)//2
        y = c**k
        if y == n: return c
        if y > n: b = c
        else: a = c
    return 0

def quadratic_check(solve_congruence):
    def wrapper(n:int, p:int):
        if n % p == 0: return {0}
        y = Legendre(n,p)
        if y == -1:
            print('NO SOLUTION : Legendre symbol (n|p)=-1 .')
            return {}
        x = powrt_int(n%p, 2)
        if x: return {x, p-x}
        return solve_congruence(n%p, p)
    return wrapper