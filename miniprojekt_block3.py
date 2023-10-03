import numpy as np


def investment(r):
    return I * ((1 + r) / r) * (((1 + r) ** n) - 1) - K


def bisection(fun, a, b, tol):
    s_fa = np.sign(fun(a))
    s_fb = np.sign(fun(b))
    assert (s_fa != s_fb)
    x = (a + b) / 2
    niter = 0
    while (b - a) / 2 > tol:
        niter = niter + 1
        s_fx = np.sign(fun(x))
        if s_fa == s_fx:
            a = x
        else:
            b = x
        x = (a + b) / 2
    return x, niter


def newton(fun, fprim, x0, tol):
    maxiter = 100
    x = x0
    dx = 2 * tol
    niter = 0
    while np.abs(dx) > tol:
        niter = niter + 1
        dx = - fun(x) / fprim(x)
        x = x + dx
        assert (niter < maxiter)
    return x, niter


def fprim(r):
    return I*(-((1+r)**n)+1+n*r*((1+r)**(-1+n))+n*(r**2)*((1+r)**(-1+n)))/(r**2)


I = 1000
n = 5
K = 6000

print("Newton:", newton(investment,fprim,0.06, 10**-6))
print("Bisection:", bisection(investment, 0.04, 0.07, 10 ** -6))
