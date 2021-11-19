import math
import matplotlib.pyplot as plt

Bucher = [[0, 0, 0, 0, 0, 0, 0], 
          [1/4, 1/4, 0, 0, 0, 0, 0], 
          [3/8, 3/32, 9/32, 0, 0, 0, 0],
          [12/13, 1932/2197, -7200/2197, 7296/2197, 0, 0, 0], 
          [1, 439/216, -8, 3680/513, -845/4104, 0, 0], 
          [1/2, -8/27, 2, -3544/2565, 1859/4104, -11/40, 0],
          [None, 16/135, 0, 6656/12825, 28531/56430, -9/50, 2/55], #5次
          [None, 25/216, 0, 1408/2565, 2197/4104, -1/5, 0] #4次
        ]
def real_f(x):
    return -math.e**(-x)/2 + 2*math.e**(-2*x)/5 + 3*math.sin(x)/10 + math.cos(x)/10

def f(x, y, z, id):
    if id == 1: return math.cos(x) - 2*y - 3*z
    else: return z

def internal(_x, _y, _dy, h, flag):
    x, y, dy = _x, _y, _dy
    F = [[0] * 6 for _ in range(2)]
    for i in range(6):
        nx = x + h * Bucher[i][0]
        for j in range(2):
            ny = y
            ndy = dy
            for k in range(1, 7):
                ny += Bucher[i][k] * F[0][k - 1]
                ndy += Bucher[i][k] * F[1][k - 1]
            F[j][i] = h*f(nx, ny, ndy, j)
    x = x + h
    for j in range(1, 7):
        dy += F[1][j - 1] * Bucher[6 + flag][j]
        y += F[0][j - 1] * Bucher[6 + flag][j]
    return x, y, dy

def Runge_Kutta_Fehlberg(h, x0, y0, dy0, EPS, alpha): 
    '''
        h     : 刻み幅
        x0    : xの初期値
        y0    : y(x0)
        dy0   : dy(x0)/dxの初期値
        EPS   : 許容誤差
        alpha : 安全係数
    '''
    rx = [x0] # record of x
    ry = [y0] # record of y
    rh = [h] # record of h
    rdy = [dy0] # record of dy/dx
    while rx[-1] < 2*math.pi:
        while True:
            x5, y5, dy5 = internal(rx[-1], ry[-1], rdy[-1], h, 0) #5次
            _, y4, _ = internal(rx[-1], ry[-1], rdy[-1], h, 1) #4次
            diff = y5 - y4
            if abs(diff) < EPS :
                rh.append(h)
                h = alpha * h * ((EPS / abs(diff))**(0.2))
                break
            h /= 2.0
        rx.append(x5)
        ry.append(y5)
        rdy.append(dy5)
    
    return rx, ry, rh

res_x, res_y, res_h = Runge_Kutta_Fehlberg(math.pi / 50, 0, 0, 0, 0.00001, 0.8)

#描画
fig, ax = plt.subplots(2, 1)
ax[0].plot(res_x, res_y, linewidth = 5, label="RungeKutta")
ax[0].plot(res_x, list(map(real_f, res_x)), label="Real-Function")
ax[0].legend()
ax[1].plot(res_x, res_h, label='step width')
ax[1].set_ylim(0, 0.21)
ax[1].legend()
plt.show()
