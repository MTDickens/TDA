import numpy as np

def reduce_matrix(matrix):
    #Returns [reduced_matrix, rank, nullity]
    if np.size(matrix)==0:
        return [matrix,0,0]
    m=matrix.shape[0]
    n=matrix.shape[1]
    def _reduce(x):
        #We recurse through the diagonal entries.
        #We move a 1 to the diagonal entry, then
        #knock out any other 1s in the same  col/row.
        #The rank is the number of nonzero pivots,
        #so when we run out of nonzero diagonal entries, we will
        #know the rank.
        nonzero=False
        #Searching for a nonzero entry then moving it to the diagonal.
        for i in range(x,m):
            for j in range(x,n):
                if matrix[i,j]==1:
                    matrix[[x,i],:]=matrix[[i,x],:]
                    matrix[:,[x,j]]=matrix[:,[j,x]]
                    nonzero=True
                    break
            if nonzero:
                break
        #Knocking out other nonzero entries.
        if nonzero:
            for i in range(x+1,m):
                if matrix[i,x]==1:
                    matrix[i,:] = np.logical_xor(matrix[x,:], matrix[i,:])
            for i in range(x+1,n):
                if matrix[x,i]==1:
                    matrix[:,i] = np.logical_xor(matrix[:,x], matrix[:,i])
            #Proceeding to next diagonal entry.
            return _reduce(x+1)
        else:
            #Run out of nonzero entries so done.
            return x
    rank=_reduce(0)
    return [matrix, rank, n-rank]

# # An alternate "reduce_matrix" implementation using sympy
# from sympy import Matrix, FF
# from sympy.matrices.normalforms import smith_normal_form
# def reduced_matrix(matrix):
#     matrix = np.matrix(smith_normal_form(Matrix(matrix), domain=FF(2)))
#     ## TODO: Get the rank
#     ## TODO: Return [matrix, rank, n-rank]

# #Initialize boundary matrices
# ## a, b, c, d -> 0
# boundaryMap0 = np.matrix([[0,0,0,0]])
# ## ab -> a + b, bc -> b + c, ca -> a + c, cd -> c + d, bd -> b + d
# boundaryMap1 = np.matrix([[1,0,1,0,0],
#                           [1,1,0,0,1],
#                           [0,1,1,1,0],
#                           [0,0,0,1,1]])
# boundaryMap2 = np.matrix([[1,1,1,0,0]])

# #Smith normal forms of the boundary matrices
# smithBM0, _, rank_z_0 = reduce_matrix(boundaryMap0)
# smithBM1, rank_b_0, rank_z_1 = reduce_matrix(boundaryMap1)
# smithBM2, rank_b_1, rank_z_2 = reduce_matrix(boundaryMap2)

# betti0 = (rank_z_0 - rank_b_0)
# betti1 = (rank_z_1 - rank_b_1)
# betti2 = 0  #There is no n+1 chain group, so the Betti is 0

# print(smithBM0)
# print(smithBM1)
# print(smithBM2)
# print("Betti #0: %s \n Betti #1: %s \n Betti #2: %s" % (betti0, betti1, betti2))


# return the n-simplices in a complex
def nSimplices(n, complex):
    nchain = [simplex for simplex in complex if len(simplex) == n+1]
    if (len(nchain) == 0): nchain = [0]
    return nchain

#check if simplex is a face of another simplex
def checkFace(face, simplex):
    if simplex == 0:
        return 1
    elif set(face) < set(simplex): #if face is a subset of simplex
        return 1
    else:
        return 0

#build boundary matrix for dimension n ---> (n-1) = p
def boundaryMatrix(nchain, pchain):
    bmatrix = np.zeros((len(nchain),len(pchain)))
    i = 0
    for nSimplex in nchain:
        j = 0
        for pSimplex in pchain:
            bmatrix[i, j] = checkFace(pSimplex, nSimplex)
            j += 1
        i += 1
    return bmatrix.T

# S = [{0}, {1}, {2}, {3}, {0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 1}, {0, 1, 2}] #this is our simplex from above

# chain2 = nSimplices(1, S)
# chain1 = nSimplices(0, S)
# print(reduce_matrix(boundaryMatrix(chain2, chain1)))

def betti(complex):
    max_dim = len(max(complex, key=len)) #get the maximum dimension of the simplicial complex, 2 in our example
    betti_array = np.zeros(max_dim) #setup array to store n-th dimensional Betti numbers
    z_n = np.zeros(max_dim) #number of cycles (from cycle group)
    b_n = np.zeros(max_dim) #b_(n-1) boundary group
    for i in range(max_dim): #loop through each dimension starting from maximum to generate boundary maps
        bm = 0 #setup n-th boundary matrix
        chain2 = nSimplices(i, complex) #n-th chain group
        if i==0: #there is no n+1 boundary matrix in this case
            bm = 0
            z_n[i] = len(chain2)
            b_n[i] = 0
        else:
            chain1 = nSimplices(i-1, complex) #(n-1)th chain group
            bm = reduce_matrix(boundaryMatrix(chain2, chain1))
            z_n[i] = bm[2]
            b_n[i] = bm[1] #b_(n-1)

    for i in range(max_dim): #Calculate betti number: Z_n - B_n
        if (i+1) < max_dim:
            betti_array[i] = z_n[i] - b_n[i+1]
        else:
            betti_array[i] = z_n[i] - 0 #if there are no higher simplices, the boundary group of this chain is 0

    return betti_array