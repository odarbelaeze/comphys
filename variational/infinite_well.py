from scipy.linalg import eig
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

def S_func(m, n):
    if ((n + m) % 2 != 0): return 0.0
    return 2.0 / (n + m + 5) - \
           4.0 / (n + m + 3) + \
           2.0 / (n + m + 1)

def H_func(m, n):
    if ((n + m) % 2 != 0): return 0.0
    return - 8.0 * (1 - m - n - 2 * m * n) / \
             ((m + n + 3) * (m + n + 1) * (m + n - 1))

def main():
    N = 10
    H = np.zeros([N, N])
    S = np.zeros([N, N])

    for m in xrange(N):
        for n in xrange(N):
            S[m][n] = S_func(m, n)
            H[m][n] = H_func(m, n)

    E, C = eig(H, S)
    print E
    ground = E.min()
    i = np.argmax(E == ground)

    x = np.arange(- 1.0, 1.0, 0.01)
    y = 0.0 * x[:]

    for p in xrange(N):
        y_p = C.T[i][p] * (x ** p) * (x - 1) * (x + 1)
        y = y + y_p
        if C.T[i][p] is not 0.0:
            plt.plot(x, y_p, '--', label = "p = %d, cp = %f" % (p, C.T[i][p]))

    y2 = np.cos(np.pi / 2.0 * x)

    plt.plot(x, y)
    plt.plot(x, y2, '-o')
    plt.legend()
    plt.grid()
    plt.show()

    print C[i], ground

if __name__ == '__main__':
    main()
