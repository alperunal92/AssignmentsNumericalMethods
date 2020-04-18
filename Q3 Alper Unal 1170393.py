from numpy import zeros,argmax,dot
from modules import swap

def gaussPivot(a,b,tol=1.0e-12):
    n = len(b)

''' ## module gaussPivot
    x = gaussPivot(a,b,tol=1.0e-9).
    Solves [a]{x} = {b} by Gauss elimination with
    scaled row pivoting
'''
  # Set up scale factors
    s = zeros(n)
    for i in range(n):
        s[i] = max(abs(a[i,:]))

    for k in range(0,n-1):

      # Row interchange, if needed
        p = argmax(abs(a[k:n,k])/s[k:n]) + k
        if abs(a[p,k]) < tol: print('Matrix is singular')
        if p != k:
            swap.swapRows(b,k,p)
            swap.swapRows(s,k,p)
            swap.swapRows(a,k,p)

      # Elimination
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a [i,k+1:n] - lam*a[k,k+1:n]
                b[i] = b[i] - lam*b[k]
    if abs(a[n-1,n-1]) < tol: print('Matrix is singular')

  # Back substitution
    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    return b

    # Main program starts here
    if __name__ == '__main__':
        a = np.array([[1,  -1,   2,  1],
                      [3,   2,   1,  4],
                      [5,   8,   6,  3],
                      [4,   2,   5,  3]])
        b = np.array([1, 1, 1, -1])
        x = gauss(a, b)
        print('Gauss result is x = \n %s' % x)
