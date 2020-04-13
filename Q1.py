from numpy 
import array, zeros, fabs, linalg

a = array([[2, -2, 0],[1, -1, 1],[0, -3, -2]], float)
b = array([-6, 1, 5], float)


def solve(a, b):
    n = len(b)
    x = zeros(n, float)
    eps = 0.5e-15

    # Elimination
    for k in range(n):
        if fabs(a[k, k]) < eps:
            for i in range(k + 1, n):
                if (fabs(a[i, k] > a[k, k])):
                    a[[k, i]] = a[[i, k]]
                    b[[k, i]] = b[[i, k]]
                    break

        for i in range(k + 1, n):
            if a[i, k] == 0:
                continue
            factor = a[k, k] / a[i, k]
            for j in range(k, n):
                a[i, j] = a[k, j] - a[i, j] * factor
            b[i] = b[k] - factor * b[i]

    # Back substitution
    x[n - 1] = b[n - 1] / a[n - 1][n - 1]

    for i in range(n - 2, -1, -1):
        sum = 0
        for j in range(i + 1, n):
            sum += (x[j] * a[i, j])
        x[i] = (b[i] - sum) / a[i, i]

    return x


def compare_with_numpy(a, b):
    numpy_sol = linalg.solve(a, b)
    solution = solve(a, b)
    print(solution)



compare_with_numpy(a, b)
