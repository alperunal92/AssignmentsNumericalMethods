import numpy as np

A = np.array([[1,  -1,   2,  1],
              [3,   2,   1,  4],
              [5,   8,   6,  3],
              [4,   2,   5,  3]])
R = np.array([[1],[1],[1],[-1]])
AR=np.append(A,R,axis=1)  #append the right-handside vector into coefficient matrix as a column
n=len(AR)      #size of the coefficient array (number of rows)

for k in range(0,n):
    pivot=abs(AR[k,k])
    p_row=k
    for i in range(k+1,n):  #find the largest element in pivot column
        if abs(AR[i,k])>pivot:
            p_row=i
            pivot=abs(AR[i,k])
    if p_row>k:
        for j in range (0,n+1):  #exchange row with biggest element with the pivot row
            temp=AR[k,j]
            AR[k,j]=AR[p_row,j]
            AR[p_row, j]=temp
    for i in range(k+1,n):  #Use pivot row to eliminate elements below main diagonal
        s=-AR[i,k]/AR[k,k]
        for j in range(k,n+1):
            AR[i,j]=AR[i,j]+s*AR[k,j]

X=np.zeros([n,1])
X[n-1]=AR[n-1,n]/AR[n-1,n-1]

for i in range(n-1,-1,-1):
    sum=0
    for j in range(i+1,n):
        sum=sum+AR[i,j]*X[j]
    X[i]=(AR[i,n]-sum)/AR[i,i]
print(X)