import utilities2 as ty
import numpy as np
import matplotlib.pyplot as plt
print("solution using L_U")
ty.L_U("L.txt","U.txt")
#ty.jacobi("L.txt","U.txt")
ty.matrix_r("L.txt")
q2 = np.array([[ 2, -3,  0,  0,  0,  0],
             [-1,  4, -1,  0, -1,  0],
             [ 0, -1,  4,  0,  0, -1],
             [ 0,  0,  0,  2, -3,  0],
             [ 0, -1,  0, -1,  4, -1],
             [ 0,  0, -1,  0, -1,  4]])


#ty.jacobi(q2,q2b,1e-4)
inverse, ite,res=ty.inv_cg(q2)
print("conguate gradient inverse",inverse)
print(res)
print(ite)
inverse_, ite_,res_=ty.jacobi_inv(q2)
print("inverse jacobi",inverse_)
plt.plot(res,ite)
plt.plot(ite_,res_)
# plt.xlabel("iterations")
# plt.ylabel("residue")
plt.show()
