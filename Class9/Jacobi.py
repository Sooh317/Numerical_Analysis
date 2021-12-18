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
    return B

res = Jacobi_algorithm(A)
print(res)
ans, ansv = np.linalg.eig(A)
print(ans)
