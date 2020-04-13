from numpy 
import array, zeros, fabs, linalg, numpy as np


def GaussElimination(M):
    """
    Gaussian elimination method for systems of equations.
    the function receives a matrix in the following way:
    f1 -->  2*x -2y   = 6
    f2 -->  x - y + z = 1
    f3 -->  3y  -2*z = -5
    M = [
        [ 2, -2,  0,   6],
        [ 1, -1,  1,   1],
        [ 0,  3, -2,  -5]
    ]
    :param M: Matrix of system of equations
    :return b: Values of the system variables of equations in a list []
    Example:
    r = GaussElimination(M)
    """

    x = len(M)
    b = []
    m = []

    for i in range(0, x):
        b.append(M[i][x])
        m.append(M[i][:x])

    n = len(b)

    for c in range(0, n - 1):
        for r in range(c + 1, n):
            if m[r][c] != 0.0:
                temp = m[r][c] / m[c][c]
                m[r][c + 1:n] -= np.asarray(temp) * m[c][c + 1:n]
                b[r] -= np.asarray(temp) * b[c]
    for c in range(n - 1, -1, -1):
        b[c] = (b[c] - np.dot(m[c][c + 1:n], b[c + 1:n])) / m[c][c]

    return b
