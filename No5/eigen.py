import numpy as np

print('eigenvalues')
for n in range(1, 11):
    a = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            a[i][j] = 1/(i + j + 1)
    values, vectors = np.linalg.eig(a)
    for i in range(len(values)): values[i] = abs(values[i])
    mx, mn = max(values), min(values)
    print('n = %d -> %f' % (n, abs(mx/mn)))
