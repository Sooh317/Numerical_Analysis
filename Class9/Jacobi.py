import numpy as np
from copy import deepcopy

A = np.array([[4.0, -2.5, -1.0, 0.5],
              [-2.5, 4.0, 0.5, -1.0],
              [-1.0, 0.5, 4.0, -2.5],
              [0.5, -1.0, -2.5, 4.0]])

def build_rotation_matrix(n, i, j, cos, sin):
    R = np.identity(n)
    R[i][i] = cos
    R[j][j] = cos
    R[i][j] = sin
    R[j][i] = -sin
    return R

def Jacobi_algorithm(A):
    N = len(A)
    B = deepcopy(A)
    V = np.identity(N)
    calc_tau = lambda A, i, j: (A[j][j] - A[i][i]) / (2*A[i][j])
    for _ in range(3):
        for i in range(N):
            for j in range(i+1,N):
                tau = calc_tau(B, i, j)
                if tau > 0: t = -1 / (tau + np.sqrt(1 + tau**2))
                else: t = -1 / (tau - np.sqrt(1 + tau**2))
                cos = 1 / np.sqrt(1 + t**2)
                sin = t / np.sqrt(1 + t**2)
                R = build_rotation_matrix(N, i, j, cos, sin)
                B = R.dot(B).dot(R.T)
                V = V.dot(R.T)
    return B, V

def check(A, res, resv):
    for i in range(4):
        print(res[i][i] * resv[:, i])
        print(A.dot(resv[:,i]))

res, resv = Jacobi_algorithm(A)
print('res\n', res)
print('resv\n', resv)
check(A, res, resv)

'''
実行結果
(ほぼ対角)行列:
[[ 1.00000000e+00 -1.41125368e-04  7.83217609e-06 -3.00805467e-08]
 [-1.41125368e-04  8.00000000e+00 -6.69492846e-08  7.21182186e-12]
 [ 7.83217609e-06 -6.69492841e-08  5.00000000e+00 -6.48892319e-20]
 [-3.00805466e-08  7.21159129e-12 -3.88894518e-16  2.00000000e+00]]
固有ベクトルの行列:
[[ 0.50001107 -0.49998991  0.49999903 -0.49999998]
 [ 0.49998896  0.50001007 -0.50000099 -0.49999998]
 [ 0.49999088  0.50001009  0.49999901  0.50000002]
 [ 0.50000909 -0.49998993 -0.50000097  0.50000002]]

check(A, res, resv)を呼び出すと、Av = λv となっていることが確認できる
つまり、対角成分が固有値になっていることが分かる
よって
固有値:[1, 8, 5, 2]
固有ベクトル:
[0.5, 0.5, 0.5, 0.5]
[-0.5, 0.5, 0.5, -0.5]
[0.5, -0.5, 0.5, -0.5]
[-0.5, -0.5, 0.5, 0.5]
'''
