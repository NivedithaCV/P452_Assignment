#import utilities as li
import utilities2 as ty
#Y,Z=ty.Matrix_read("A.txt","B.txt")
#A,B=ty.Matrix_read("A.txt","B.txt")
#C,D = li.GaussJordan(A, B)
R,M=ty.G_J("A.txt","B.txt")
print("Using Guass_Jordan")
print("solution is",M)

print("Using L_U")
ty.L_U("A.txt","B.txt")

#output Using Guass_Jordan
# solution is [[-1.7618170439978567], [0.8962280338740136], [4.051931404116157], [-1.6171308025395428], [2.041913538501914], [0.15183248715593495]]
# Using L_U
# matrix Y [[19.0], [2.0], [-6.399999999999999], [-214.42857142857142], [-12.25], [-7.842805005213755]]
# solution matrix is [[-1.761817043997862], [0.8962280338740133], [4.051931404116158], [-1.6171308025395421], [2.041913538501913], [0.15183248715593525]]
