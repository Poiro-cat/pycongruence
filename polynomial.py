import random
from .utils import extend_gcd

class Polynomial():
    def __init__(self, coef={1:1}, modulo:int=0):
        if type(coef) == int: self.coef = {0:coef} if coef != 0 else {}
        elif type(coef) == list:
            self.coef = {}
            for i,c in enumerate(coef):
                if c != 0: self.coef[len(coef)-1-i] = c
        elif type(coef) == dict: self.coef = coef
        else: raise TypeError('Input coef must be int, list or dictionary.')
        self.degree = max(self.coef) if self.coef else -1
        self.modulo = modulo
    def simplify(self):
        for k in list(self.coef):
            self.coef[k] %= self.modulo
            if not self.coef[k]: del self.coef[k]
        self.degree = max(self.coef) if self.coef else -1
    def copy(self):
        return Polynomial(self.coef.copy(), self.modulo)
    def __add__(self, other):
        res = self.coef.copy()
        if type(other) == int:
            if 0 in res: res[0] += other
            else: res[0] = other
            p = Polynomial(res,self.modulo)
        elif type(other) == Polynomial:
            assert self.modulo == other.modulo, 'Modulo numbers are different.'
            c = other.coef
            for k in c:
                if k in res: res[k] += c[k]
                else: res[k] = c[k]
            p = Polynomial(res, self.modulo)
        else: raise TypeError('Operated instances must be Integer or Polynomial.')
        p.simplify()
        return p
    def __radd__(self, a:int):
        return self + a
    def __neg__(self):
        res = self.coef.copy()
        for k in res: res[k] = -res[k]
        return Polynomial(res, self.modulo)
    def __sub__(self, other):
        return self + (-other)
    def __rsub__(self, other):
        return -self + other
    def __mul__(self, other):
        if type(other) == int:
            res = self.coef.copy()
            for k in res: res[k] *= other
            p = Polynomial(res, self.modulo)
        elif type(other) == Polynomial:
            assert self.modulo == other.modulo, 'Modulo numbers are different.'
            c1,c2,res = self.coef, other.coef, {}
            for k1 in c1:
                for k2 in c2:
                    if k1+k2 in res: res[k1+k2] += c1[k1]*c2[k2]
                    else: res[k1+k2] = c1[k1]*c2[k2]
            p = Polynomial(res, self.modulo)
        else: TypeError('Operated instances must be Integer or Polynomial.')
        p.simplify()
        return p
    def __rmul__(self, a:int):
        return self * a
    def __pow__(self, n:int):
        res = Polynomial(1, self.modulo)
        p = self.copy()
        while True:
            if n % 2 > 0: res *= p
            n //= 2
            if n == 0: break
            p *= p
        return res
    def divmod(self, other):
        if type(other) != Polynomial: raise TypeError('Operated instances must be Polynomial.')
        assert self.modulo == other.modulo, 'Modulo numbers are different.'
        assert other.coef, 'Modulo polynomial cannot be ZERO.'
        if not self.coef: return self,self
        n,m = self.degree, other.degree
        quotient = {}
        remainder = self.copy()
        c = other.coef
        while n >= m:
            r = extend_gcd(remainder.coef[n], c[m])[0]
            r1, r2 = remainder.coef[n]//r, c[m]//r
            r, r2, _ = extend_gcd(r2, self.modulo)
            if r > 1: raise ValueError(r)           # nontrivial factor of N is found
            quotient[n-m] = r1*r2
            remainder -= other*r1*r2 * Polynomial({n-m:1}, self.modulo)
            n = remainder.degree
        quotient = Polynomial(quotient, self.modulo)
        quotient.simplify()
        return quotient, remainder
    def __mod__(self, other):
        return self.divmod(other)[1]
    def __floordiv__(self, other):
        return self.divmod(other)[0]
    def __str__(self):
        if not self.coef: return '0'
        pows = list(self.coef)
        pows.sort(reverse=True)
        st = ''
        for k in pows:
            c = self.coef[k]
            if not c: continue
            if c > 0: st = st + ' + '
            if c < 0: st = st + ' - '
            if not k or c*c > 1: st = st + str(c if c > 0 else -c)
            if c*c > 1 and k: st = st + '*'
            if k: st = st + 'X'
            if k > 1: st = st + f'^{k}'
        st = st[3:] if st[1]=='+' else st[1:]
        return st
    def monic(self):
        assert self.modulo > 0
        if self.coef[self.degree] == 1: return self
        p = self * extend_gcd(self.coef[self.degree], self.modulo)[1]
        if p.coef[p.degree] > 1: raise ValueError(p.coef[p.degree])
        return p
    def simplify_by_fermat(self):
        assert self.modulo > 0
        for d in list(self.coef):
            if d >= self.modulo:
                j = d % (self.modulo-1)
                if j in self.coef: self.coef[j] += self.coef[d]
                else: self.coef[j] = self.coef[d]
                del self.coef[d]
        self.simplify()
    def differentition(self):
        d_coef = {}
        for i in self.coef:
            c = self.coef[i]*i % self.modulo
            if c: d_coef[i-1] = c
        return Polynomial(d_coef, self.modulo)
    def value(self, a:int):
        p = Polynomial({0:a}, self.modulo)
        p = poly_nest(self, p)
        if 0 in p.coef: return p.coef[0]
        return 0

def poly_pow_mod(p1:Polynomial, N:int, p2:Polynomial):                # computing p1(x)^N mod p2(x)
    assert p1.modulo == p2.modulo, 'Modulo numbers are different.'
    assert p2.coef, 'Modulo polynomial cannot be ZERO.'
    res = Polynomial(1,p1.modulo)
    p1 %= p2
    while True:
        if N % 2 > 0: res = res*p1 % p2
        N //= 2
        if N == 0: break
        p1 = p1*p1 % p2
    return res

def gcd_poly(p1:Polynomial, p2:Polynomial):                                 # computing greatest common divisor of two polynomials
    if not (p1.coef and p2.coef): return p1 + p2
    p1,p2 = p2, p1%p2
    while p2.coef: p1,p2 = p2, p1%p2
    return p1

def poly_nest(p1:Polynomial, p2:Polynomial):                  # computing p1(p2(x))
    assert p1.modulo == p2.modulo, 'Modulo numbers are different.'
    p3 = Polynomial({}, p1.modulo)
    y = Polynomial({0:1}, p1.modulo)
    for i in range(p1.degree+1):
        if i in p1.coef: p3 += p1.coef[i]*y
        y *= p2
    return p3
