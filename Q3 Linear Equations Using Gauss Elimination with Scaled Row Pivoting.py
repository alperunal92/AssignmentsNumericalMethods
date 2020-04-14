import numpy as np

def gauss_pivot(a,b,tol=1.0e-12):
    """
    x = gaussPivot(a,b,tol=1.0e-12).
    Solves [a]{x} = {b} by Gauss elimination with
    scaled row pivoting
    """
    a = np.copy(a)
    b = np.copy(b)
    n = len(b)
    assert (np.all(np.shape(a) ==(n,n))) # check if a is a square matrix
    # Set up scale factors
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(a[i,:])) # find the max of each row
    for k in range(0, n-1): #pivot row
        # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) # find which row has max item for each col k, and scale by s
        if abs(a[p,k]) < tol: 
            raise Exception("Matrix is singular")
        if p != k: # swap rows if current row does not contain max item with the one contains max item within same col
            a[[k,p+k],:] = a[[p+k, k],:]
            b[k],b[p+k] = b[p+k],b[k] 
            s[k],s[p+k] = s[p+k],s[k]

        # Elimination phase of matrix a
        for i in range(k+1,n):
            if a[i,k] != 0.0: # skip if a(i,k) is already zero
                lam = a [i,k]/a[k,k] 
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                b[i] = b[i] - lam*b[k]

    if abs(a[n-1,n-1]) < tol: 
        raise Exception("Matrix is singular")

    # Back substitution phase, solution is substituted by b

    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
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