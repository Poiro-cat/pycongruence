import random
from .utils import (
    extend_gcd, power_mod, power_mod_quadratic_field, Legendre, trial_div_factorize, Sunzi, quadratic_check
)
from .polynomial import Polynomial, poly_pow_mod, gcd_poly, poly_nest


##### 1. Linear equation #####

def Solve_linear(a:int, b:int, N:int):
    assert a % N > 0, 'Modulus N should not divide input a.'
    r, m, _ = extend_gcd(a,N)
    if b % r > 0:
        #print('NO SOLUTION : Input b can not be divided by gcd(a,N).')
        return set()
    return {(m*b//r)%N}


##### 2. Quadratic equation mod Prime #####

@ quadratic_check
def Solve_quadratic_mod_prime_Cipolla(n:int, p:int):
    while True:
        a = random.randint(0, p-1)
        D = (a*a - n) % p
        if D == 0: return {a, p-a}
        if Legendre(D,p) == 1: continue
        x, y = power_mod_quadratic_field(a, 1, D, (p+1)//2, p)
        assert (x*x-n)%p == y == 0, 'Modulus p is not prime (Fermat\'s little theorem failed).'
        return {x, p-x}

@ quadratic_check
def Solve_quadratic_mod_prime_Tonelli_Shanks(n:int, p:int):
    r, s = 0, p-1
    while s % 2 == 0:
        r += 1
        s //= 2
    while True:
        z = random.randint(1, p-1)
        if Legendre(z,p) == -1: break
    #print('z =', z)
    X, S, T, rho = power_mod(n, (s+1)//2, p), power_mod(n,s,p), power_mod(z,s,p), r-1
    while True:
        #print(f'X = {X} ,  S = {S},  T = {T} ,  rho = {rho}')
        if S == 1: return {X, p-X}
        rho_, S_ = 0, S
        while S_ != 1:
            S_ = S_*S_ % p
            rho_ += 1
        for _ in range(rho-rho_): T = T*T % p
        rho = rho_
        X = X*T % p
        S = S*T*T % p

@ quadratic_check
def Solve_quadratic_mod_prime_Pocklington(n:int, p:int):
    if p % 4 == 3:
        x = power_mod(n, p//4+1, p)
    elif p % 8 == 5:
        k = p // 8
        x = power_mod(n, k+1, p)
        if power_mod(n, 2*k+1, p) == p-1: x = x * power_mod(2, 2*k+1, p) % p
    else:
        while True:
            t, u = random.randint(0, p-1), random.randint(0, p-1)
            W = (t*t + n*u*u) % p
            if Legendre(W,p) == -1: break
        r, s = 0, p-1
        while s % 2 == 0:
            r += 1
            s //= 2
        tm, um = power_mod_quadratic_field(t,u,-n,s,p)
        assert tm*um > 0, 'Sequence ts=0 or us=0 : contradictory to that W^s is a quadratic non-residue.'
        for i in range(1, r+1):
            if (tm*tm - n*um*um) % p == 0:
                x = list(Solve_linear_congruence(um,tm,p))[0]
                break
            assert i < r, 'Sequence t_{2^i*s}=0 not found for 1<=i<r.'
            tm, um = (tm*tm-n*um*um) % p, 2*tm*um % p
    return {x, p-x}


##### 3. Polynomial equation mod Prime #####

def Solve_polynomial_mod_prime(f, p:int):
    if type(f) in (list,dict): f = Polynomial(f)
    assert type(f) == Polynomial
    f.modulo = p
    f.simplify_by_fermat()
    #print(f'\nInput f(x) = {f}')
    if f.degree == -1: return {x for x in range(p)}
    dmin = min(f.coef)
    solutions = set() if dmin == 0 else {0}
    x = Polynomial(modulo=p)
    if dmin > 0: f //= x**dmin
    if f.degree == 0: return solutions
    if f.degree == 1:
        a, b = f.coef[1], -f.coef[0]
        return solutions | Solve_linear(a,b,p)
    F0 = poly_pow_mod(x, (p-1)//2, f)
    #print(f'factors:   {gcd_poly(f,F0-1)}  |  {gcd_poly(f,F0+1)}')
    for F in [F0-1, F0+1]:
        f0 = gcd_poly(f,F)
        z = random.randint(1, p-1)
        ys = Solve_polynomial_mod_prime(poly_nest(f0, x-z), p)
        solutions |= {(y-z)%p for y in ys}
    return solutions


##### 4. Polynomial equation mod Prime Power #####

def Solve_polynomial_mod_prime_power(f, p:int, m:int):
    if type(f) in (list,dict): f = Polynomial(f)
    assert type(f) == Polynomial
    f.modulo = p**m
    f.simplify()
    S = Solve_polynomial_mod_prime(f.copy(),p)
    #print(f'\n##  j = 1  |  S = {S}\n')
    df = f.copy()
    df.modulo = p
    df = df.differentition()
    pj = p
    for j in range(2, m+1):
        if not S: return S
        S_ = set()
        pj *= p
        for xi in S:
            dfxi = df.value(xi)
            fxi = f.value(xi) % pj
            #print(f'    xi = {xi}, fxi = {fxi}, dfxi = {dfxi}')
            if dfxi != 0:
                lmd = fxi // (pj//p)
                t = list(Solve_linear(dfxi, -lmd, p))[0]
                S_.add(xi + t * (pj//p))
            elif fxi == 0:
                for t in range(p): S_.add(xi + t * (pj//p))
        S = S_.copy()
        #print(f'\n##  j = {j}  |  S = {S}\n')
    return S


##### 5. General case + Equation group #####

def Solve_congruence(*args):
    eqs = {}
    for f,N in args:
        factorization = trial_div_factorize(N)
        for p in factorization:
            m = factorization[p]
            print(f'\n###  p**m = {p}**{m}  |  f = {f}')
            ans = Solve_polynomial_mod_prime_power(f,p,m)
            print('     ans =', ans)
            if not ans: return ans
            if p not in eqs: eqs[p] = [m, ans]
            else:
                pd = p**abs(eqs[p][0]-m)
                p0 = p**min(eqs[p][0], m)
                a0, a1 = (eqs[p][1], ans) if m > eqs[p][0] else (ans, eqs[p][1])
                eqs[p][0] = max(eqs[p][0], m)
                eqs[p][1] = a1 & {xi+t*p0 for xi in list(a0) for t in range(pd)}
            print(eqs)
    Q = 1
    solutions = set()
    for p in eqs:
        if not eqs[p][1]: return set()
        q = p**eqs[p][0]
        if Q == 1: solutions = eqs[p][1].copy()
        else: solutions = {Sunzi([a, b], [Q, q])[0] for a in solutions for b in eqs[p][1]}
        Q *= q
        print(f'p**m = {q} ,  ans = {eqs[p][1]}  |  sol = {solutions}')
    return solutions, Q