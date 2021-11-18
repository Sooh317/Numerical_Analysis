import math

def f(x): # 目的の関数
    return (math.e**(x**4)) - math.pi

def e(sign, x, p): #補助関数(e^(sign * x^p))
    return math.e**(sign * (x**p))

def c(x): # f/f' = (1 - pi * e^(-x^4)) / 4x^3
    return (1.0 - math.pi * e(-1, x, 4)) / (4*(x**3))

def g(x, lamb, j):
    r = lambda x, lamb, j: (x - lamb**(-j)/(4*(x**3)))**4 - x**4
    y = lambda x, lamb, j: (x - lamb**(-j)*c(x))
    epsilon = 10**(-18)
    return r(x, lamb, j) + math.log(abs(1 - math.pi * e(-1, y(x, lamb, j), 4)) + epsilon) - math.log(abs(1 - math.pi * e(-1, x, 4)) + epsilon)



def Newton(a): #f(x) = 0 を解く, 初期値をaでとく
    rx = [a] # record of x
    EPS = 10**(-9)
    while True:
        nx = rx[-1] - c(rx[-1])
        error = abs(nx - rx[-1])
        rx.append(nx)
        if error < EPS : break
    return rx

def decelaration_Newton(a, beta, lamb):
    rx = [a] # record of x
    EPS = 10**(-9)
    while True:
        j = -1
        nx = -1.0
        for i in range(0, 10000000):
            nx = rx[-1] - (lamb**(-i))*c(rx[-1])
            L = abs(g(rx[-1], lamb, i))
            R = (1 - beta*(lamb**(-i)))
            if L <= R:
                j = i
                break
        error = abs(nx - rx[-1])
        rx.append(nx)
        if error < EPS : break
    return rx


NX_middle = Newton(0.5)
NX_small = Newton(0.0001)
NX_large = Newton(5)
print('NX_small\n', NX_small, len(NX_small))
print('NX_middle\n', NX_middle, len(NX_middle))
print('NX_large\n', NX_large, len(NX_large))
DNX_small = decelaration_Newton(0.0001, 0.5, 1.5)
DNX_middle = decelaration_Newton(0.5, 0.5, 1.5)
DNX_large = decelaration_Newton(5, 0.5, 1.5)
print('DNX_small\n', DNX_small, len(DNX_small))
print('DNX_middle\n', DNX_middle, len(DNX_middle))
print('DNX_large\n', DNX_large, len(DNX_large))
