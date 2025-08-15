# pycongruence

This is a Python 3 package providing methods for solving congruence equiations of the form $f(x)\equiv0\pmod N$, where $f(x)$ is a polynomial.

## List of methods

### 1. Solving Linear Equations

To solve linear congruences of the form $ax\equiv b\pmod N$ where $a\neq 0\pmod N$, the package includes the method:

- `Solve_linear(a, b, N)` 

### 2. Solving Quadratic Equation with Prime Moduli

For quadratic congruences of the form $x^2\equiv n\pmod p$ where $p$ is a prime number, the following three methods are provided:

- `Solve_quadratic_mod_prime_Cipolla(n, p)`: Implements Cipolla's algorithm, a straightforward approach for solving quadratic congruences.
- `Solve_quadratic_mod_prime_Tonelli_Shanks(n, p)` : Implements Tonelli_Shanks algorithm, which is based on Euler's criterion.
- `Solve_quadratic_mod_prime_Pocklington(n, p)`: Implements Pocklington's algorithm.

### 3. Solving Polynomial Equation with Prime Modulus

To solve polynomial congruences with the form $f(x)\equiv0\pmod p$ where $p$ is a prime number and $f(x)$ is an arbitrary polynomial, the package offers:
- `Solve_polynomial_mod_prime(f, p)` : Uses the Berlekamp-Rabin algorithm. The input `f` could be a list of coefficients, a dictionary of coefficients, or a `Polynomial` instance defined within `pycongruence`.


### 4. Solving Polynomial Equation with Prime Power Modulus

For polynomial congruences of the form $f(x)\equiv0\pmod{p^m}$ where $m> 1$, $p$ is a prime number, and $f(x)$ is an arbitrary polynomial, the following method is provided:
- `Solve_polynomial_mod_prime_power(f, p, m)`: The input `f` could be a list of coefficients, a dictionary of coefficients, or a `Polynomial` instance defined within `pycongruence`.

### 5. Solving Any Congruence Equations

The function `Solve_congruence(*args)`  handles arbitrary polynomial congruences of the form $f(x)\equiv0\pmod N$ where $N>1$. It also supports solving systems of congruence equations.

The input `*args` should be formatted as `[f1, N1], [f2, N2], ...`, where `f1,f2,...` could be lists of coefficients, dictionaries of coefficients, or `Polynomial` instances defined in `pycongruence`.
