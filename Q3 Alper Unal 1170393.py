import numpy as np

def zamienW(M, i,j):
    
    if len(M.shape) == 1:
        M[i],M[j] = M[j],M[i]
    else:
        M[[i,j],:] = M[[j,i], :]

def Gauss_scaled(A, B):
    
    n = len(A)
    S = np.zeros(n)         
    P = np.array(range(n))  
    A= A.astype(np.float64)
    B= B.astype(np.float64)

    
    for i in range(n):
        S[i]= max(abs(A[i,:]))

    for k in range(0, n-1):
       
        for j in range(k+1, n):
            if np.abs(A[j][k])/S[j] >= abs(A[i][k])/S[i]:
                pivot =j

        if pivot != k:
            zamienW(B,k,pivot)
            zamienW(S,k,pivot)
            zamienW(A,k,pivot)

        for i in range(k+1, n):
            A[k][k] = 0.0 if np.isclose(A[k][k], 0) else A[k][k]
            z = A[i][k] / A[k][k]
            A[i][k] = z

            for j in range(k+1, n):
                A[i][j] = A[i][j] - z* A[k][j]
            B[i] = B[i] - z*B[k]

    B[n-1] = B[n-1]/A[n-1,n-1]
    for k in range(n-2,-1,-1):
        A[k][k] = 0.0 if np.isclose(A[k][k], 0) else A[k][k]
        B[k] = (B[k] - np.dot(A[k,k+1:n],B[k+1:n]))/A[k,k]

    return B

    # Main program starts here
    if __name__ == '__main__':
        A = np.array([[1,  -1,   2,  1],
                      [3,   2,   1,  4],
                      [5,   8,   6,  3],
                      [4,   2,   5,  3]])
        B = np.array([1, 1, 1, -1])
        x = Gauss_scaled(A, B)
        print('Gauss result is x = \n %s' % x)
