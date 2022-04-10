import numpy as np
from fraction import Fraction
#A=0,2,5           B=1
#  3,-1,2           -2
#  1,-1,3            3
def G_J(a_file,b_file):
    fhand=open(a_file)
    ghand=open(b_file)

    M=[]
    N=[]
    A=[[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    B=[[0],[0],[0],[0],[0],[0]]
#A=0,2,5            B=1
#  3,-1,2             -2
#  1,-1,3              3
    for line in fhand:
        line=line.rstrip()
        li=list(line.split(","))
        c=len(li)
        M.append(li)
    for k in ghand:
        k=k.rstrip()
        lis=list(k.split(","))
        d=len(lis)
        N.append(lis)
    r=len(M)

    for i in range(r):
        B[i][0]=int(N[i][0])
        for j in range(r):
            A[i][j]=int(M[i][j])

    def partial_pivot(a,b):
        r=len(a)
        for i in range(r):
            if a[i][i]==0:
                for k in {0,1,2}:
                    if k==i or a[i][i]!=0:
                        continue
                    else:
                        if abs(a[k][i])>abs(a[i][i]):
                            c=b[i][0]
                            b[i][0]=b[k][0]
                            b[k][0]=c
                            for j in range(r):
                                pivot=a[i][j]
                                a[i][j]=a[k][j]
                                a[k][j]=pivot

        return a,b

#Gauss Jordan method
    def Gauss_Jordan(a,b):
        """Gauss Jordan method of decomposition"""
        for q in range(r):
            pivot=a[q][q]
            b[q][0]=b[q][0]/pivot
            for l in range(q,r):
                a[q][l]= a[q][l]/pivot

            for w in range(r):
                if a[w][q]==0 or q==w:
                    continue
                else:
                    factor=a[w][q]
                    b[w][0]=b[w][0]-factor*b[q][0]
                    for c in range(q,r):
                        a[w][c]=a[w][c]-factor*a[q][c]

        return a,b

    X,Y=partial_pivot(A,B)

    A,B=Gauss_Jordan(X,Y)
    return (A,B)
#________________________________________________________________________________

def matrix_r(A):
    with open(A,"r") as fhand:
        M=[]
        N=[]
        for line in fhand:
          line=line.rstrip()
          li=list(line.split(","))
          c=len(li)
          M.append(li)
        r=len(M)
        c=len(M[0])
        A=[[0 for y in range(c)]for x in range(r)]

        for i in range(r):
          for j in range(c):
              A[i][j]=float(M[i][j])
        return(A,c)


def partial_pivot(a,b,col):
    """This function does partial pivoting of passed matrices"""
    r=len(a)
    for i in range(r):
        if a[i][i]==0:
            for k in range(i,col):
                if k==i or a[i][i]!=0:
                    continue
                else:
                    if abs(a[k][i])>abs(a[i][i]):
                        # c=b[i][col-1]
                        # b[i][col-1]=b[k][col-1]
                        # b[k][col-1]=c
                        for j in range(r):
                            pivot=a[i][j]
                            a[i][j]=a[k][j]
                            a[k][j]=pivot
                        for z in range(col):
                            c=b[i][z]
                            b[i][z]=b[k][z]
                            b[k][z]=c
    return a,b
def L_U(a_file,b_file):
    A,c_A=matrix_r(a_file)
    B,c_B=matrix_r(b_file)
    A,B=partial_pivot(A,B,c_B)
    # A=1,0,1,2         B=6
    #   0,1,-2,0         -3
    #   1,2,-1,0         -2
    #   2,1,3,-2         0


    def L_Udec(A):
        for j in range(c_A):
            for i in range(len(A)):

                #diagonal
                if i==j:
                    sum=0
                    for u in range(i):
                        sum=sum+A[i][u]*A[u][i]
                    A[i][i]=A[i][i]-sum

                    #elements of upper triangle
                if i<j:
                    sum=0
                    for k in range(i):
                        sum=sum+A[i][k]*A[k][j]
                    A[i][j]=A[i][j]-sum

                    #elements of lower triangle
                if i>j:
                    sum=0
                    for z in range(j):
                        sum=sum+A[i][z]*A[z][j]
                    A[i][j]=(A[i][j]-sum)/A[j][j]
        return(A)

    def forw_backw(A,B,col):
        for k in range(col):
            r=len(A)

            #forward substitution
            Y=[[0] for y in range(r)]
            for i in range(r):
                sum=0
                for k in range(i):
                    sum=sum+A[i][k]*Y[k][0]
                Y[i][0]=B[i][0]-sum
            print("matrix Y",Y)

            #backward substitution
            X=[[0] for w in range(r)]
            for l in range(r-1,-1,-1):
                sum=0
                for m in range(l+1,r):
                    sum=sum+A[l][m]*X[m][0]
                X[l][0]=(Y[l][0]-sum)/A[l][l]
        print("solution matrix is",X)

    L_Udec(A)
    forw_backw(A,B,c_B)
#
def resi(x0, x_res):
    sum = 0
    for i in range(len(x0)):
        for j in range(len(x0[i])):
            x_x = abs(x_res[i][j] - x0[i][j])
            sum = sum + x_x
    return sum
# def jacobi(A,B,e,b=[[1],[2],[3]]):
#     n=len(b)
#     b_ =[[0] for x in range(len(A))]
#     def resi(b,b_):
#         sum =0
#         for i in range(len(b)):
#             for j in range(len(b[i])):
#                 difference= abs(b_[i][j]-b[i][j])
#                 sum = sum+difference
#         return(sum)
#     cycle=0
#     while resi(b,b_)>=e:
#         if cycle != 0:
#             for i in range(len(b_)):
#                 for j in range(len(b_[i])):
#                     b[i][j] = b_[i][j]
#             for i in range(len(A)):
#                 sum = 0
#                 for j in range(len(A[i])):
#                     if j != i:
#                         sum = sum + (A[i][j] * b[j][0])
#                 b_[i][0] = (1/A[i][i]) * (B[i][0] - sum)
#         cycle = cycle + 1
#     return {'X':x_res, 'iterations':count}
# def Jacobi_inverse(A,b=None):
#     r= len(A)
#     R=[[0 for i in range(r)] for j in range(r)]
#     for k in range(r):
#         C=[[0] for i in range(r)]
#         C[k][0]=1
#         D=jacobi(A,C,0.1,b)['X']
#         for l in range(r):
#             R[i][j]= round(D[j][0],2)
#
#         return(R)

def jacobi(A, B, eps, x0=[[0],[0],[0]]):

	# |x_k - x_k+1|

	iteration=[]
	res=[]
	count = 0
	x_res = [[0] for x in range(len(A))]
	while resi(x_res,x0) >= eps:
		if count != 0:
			for i in range(len(x_res)):
				for j in range(len(x_res[i])):
					x0[i][j] = x_res[i][j]
		for i in range(len(A)):
			sum = 0
			for j in range(len(A[i])):
				if j != i:
					sum = sum + (float(A[i][j]) * float(x0[j][0]))
			x_res[i][0] = (1/A[i][i]) * (B[i][0] - sum)
		iteration.append(count)
		res.append(resi(x_res,x0))
		count = count + 1

	return x_res,iteration,res

def jacobi_inv(A ):
    x0=[[2],[3],[1],[2],[1],[2]]
    A=np.asarray(A)
    n = len(A)

    inverse= [[0 for x in range(n)]for y in range(n)]

    for i in range(n):
    	b = [[0] for x in range(n)]
    	b[i][0] = 1

    	col,iteration,res = jacobi(A,b,0.1,x0)
    	for j in range(n):
    		inverse[i][j] = round(col[j][0],2)

    return inverse,iteration,res

def gauss_sidel(A,B):
    N = len(A)
    x = [0 for x in range(len(B))]

    for i in range(1, N+1):
        x_new = [0 for x in range(len(B))]
        for i in range(A.shape[0]):
            s1 = np.dot(A[i, :i], x_new[:i])
            s2 = np.dot(A[i, i + 1:], x[i + 1:])
            x_new[i] = (B[i] - s1 - s2) / A[i, i]
        if np.linalg.norm(np.subtract(x,x_new)) < 1e-5:
            break
        x = x_new
    return x

def guass_s_inverse(A,tol):

    B = [[0] for j in range(len(A))]

    matrix = np.zeros((len(A), len(A)))  # Pre-allocate matrix
    for i in range(len(A)):
        B[i]==1
        column = gauss_sidel(A, B)
        for j in range(1, len(A)):
            matrix[:, i] = column[0]
    print(column)
    return ings
def inv_g_s(A):
    a=A
    aj = np.array(A)
    bij=np.identity(len(A))
    marix = np.zeros((len(a), len(a)))  # Pre-allocate matrix
    for i in range(6):
        column,error,iteration = gauss_seidel(aj, bij[i], 1e-4)
        for j in range(1, 6):
            matrix[:, i] = column[0]
    return(matrix,error,iteration)


def conjugate_gradient(A,B,x):
    # residual
    r = B - np.matmul(A,x)
    p = r
    error=[]
    iteration=[]
    rs_old = np.matmul(np.transpose(r),r)
    for i in range(len(B)):
        Ap = np.matmul(A,p)
        alpha = rs_old/(np.matmul(np.transpose(p),Ap))
        x = x + alpha*p
        r = r - alpha*Ap
        rs_new = np.matmul(np.transpose(r),r)
        iteration.append(i)
        error.append(np.sqrt(np.sum((rs_new))))
        if np.sqrt(rs_new) < 1e-5:
            break
        else:
            p = r + (rs_new/rs_old)* p
            rs_old = rs_new
    return x,error,iteration


def inv_cg(A):
    error=[]
    iteration=[]

    A=np.asarray(A)
    n = len(A)

    inverse= [[0 for x in range(n)]for y in range(n)]
    for i in range(n):
    	b = [[0] for x in range(n)]
    	b[i][0] = 1

    	col,iteration,res = conjugate_gradient(A,b,[[1],[1],[1],[1],[1],[1]])
    	for j in range(n):
    		inverse[i][j] = round(col[j][0],2)

    return inverse,iteration,res
