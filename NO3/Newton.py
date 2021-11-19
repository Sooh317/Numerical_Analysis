import math
import matplotlib.pyplot as plt

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

NX_small = Newton(0.0001)
NX_middle = Newton(0.5)
NX_large = Newton(5)
DNX_small = decelaration_Newton(0.0001, 0.5, 1.5)
DNX_middle = decelaration_Newton(0.5, 0.5, 1.5)
DNX_large = decelaration_Newton(5, 0.5, 1.5)

#描画
xs = [0.01 * i for i in range(0, 200)]
ys = list(map(f, xs))

xl, xr = 0, 6
yl, yr = -5.0, 100

fig, ax = plt.subplots(1, 3)

for i in range(3):
    ax[i].set_xlim(xl, xr)
    ax[i].set_ylim(yl, yr)
    ax[i].plot([0.01 * i for i in range(0, int(xr / 0.01))], [0 for i in range(0, int(xr / 0.01))], c = 'black')
    ax[i].plot(xs, ys)
    ax[i].set_ylabel('y')
    ax[i].set_xlabel('x')

# small
ax[0].set_title('small')
ax[0].scatter(NX_small, [0 for _ in range(len(NX_small))], marker = "o", s = 20, c = 'red', label='newton')
ax[0].scatter(DNX_small, [0 for _ in range(len(DNX_small))], marker = "^", s = 10, c = 'black', label = 'D-newton')
#middle
ax[1].set_title('middle')
ax[1].scatter(NX_middle, [0 for _ in range(len(NX_middle))], marker = "o", s = 20, c = 'red', label='newton')
ax[1].scatter(DNX_middle, [0 for _ in range(len(DNX_middle))], marker = "^", s = 10, c = 'black', label = 'D-newton')
#large
ax[2].set_title('large')
ax[2].scatter(NX_large, [0 for _ in range(len(NX_large))], marker = "o", s = 20, c = 'red', label='newton')
ax[2].scatter(DNX_large, [0 for _ in range(len(DNX_large))], marker = "^", s = 10, c = 'black', label = 'D-newton')

plt.legend()
plt.show()


print('NX_small:', len(NX_small))
print('NX_middle:', len(NX_middle))
print('NX_large:', len(NX_large))
print('DNX_small:', len(DNX_small))
print('DNX_middle:', len(DNX_middle))
print('DNX_large:', len(DNX_large))
