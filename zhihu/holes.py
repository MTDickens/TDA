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

n = 30 # Number of points to generate

# Generate space of parameter
theta = np.linspace(0, 2.0 * np.pi, n)

# Parameters of the circle
# (a, b) is the center of the circle
# r is the radius of the circle
a, b, r, = 0.0, 0.0, 5.0

x = a + r * np.cos(theta)
y = b + r * np.sin(theta)

# # Visualization of the circle
# plt.plot(x, y)
# plt.show()

# Add gaussian noise to the points
mean, std_div = 0, .35
x2 = np.random.normal(mean, std_div, n) + x
y2 = np.random.normal(mean, std_div, n) + y

# plt.scatter(x_gauss, y_gauss)
# plt.show()

newData = np.array(list(zip(x2,y2)))
import SimplicialComplex

eps = 6
graph = SimplicialComplex.buildGraph(raw_data=newData, epsilon=eps)
ripsComplex = SimplicialComplex.rips(graph=graph, k=3)

# Genus = 1
# Still bugs here :(
lower = [{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}] + [{0,5},{1,6},{2,7},{3,4}]
upper = [set(map(lambda x: x+8, a)) for a in lower]
connect = [{i, i+8} for i in range(8)] + [{i, (i+9) % 16} for i in range(8)]
lower_faces = [{0,1,5},{0,4,5},{1,5,6},{1,6,2},{2,3,7},{2,6,7},{3,7,4},{0,3,4}]
upper_faces = [set(map(lambda x: x+1, a)) for a in lower_faces]
connect_upper_triangle = [{i % 16, (i+1) % 16, (i+9) % 16} for i in range(8)]
connect_lower_triangle = [{i % 16, (i+8) % 16, (i+9) % 16} for i in range(8)]
debug = []
# debug = [{4,5,6},{4,6,7}, {4+8,5+8,6+8},{4+8,6+8,7+8}]
ripsComplex = [{i} for i in range(16)] + lower + upper + connect + connect_upper_triangle + connect_lower_triangle + debug

betti_list = betti(ripsComplex)
print(betti_list)
V = 16
E = len(lower + upper + connect)
F = len(lower_faces + upper_faces + connect_lower_triangle + connect_upper_triangle + debug)
print("V, E, F = {}, {}, {}".format(V, E, F))
print(f"V - E + F = {V - E + F}")
print(f"b0 - b1 + b2 = {betti_list[0] - betti_list[1] + betti_list[2]}")

SimplicialComplex.drawComplex(origData=newData, ripsComplex=ripsComplex)
