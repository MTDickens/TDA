#This set of functions allows for building a Vietoris-Rips simplicial complex from point data

import numpy as np
import matplotlib.pyplot as plt

def euclidianDist(a,b):
  return np.linalg.norm(a - b) #euclidian distance metric

#Build neighorbood graph
def buildGraph(raw_data, epsilon = 3.1, metric=euclidianDist): #raw_data is a numpy array
  data_len = raw_data.shape[0]
  nodes = [x for x in range(data_len)] #initialize node set, reference indices from original data array
  edges = [] #initialize empty edge array
  weights = [] #initialize weight array, stores the weight (which in this case is the distance) for each edge
  for i in range(data_len): #iterate through each data point
      for j in range(data_len-i): #inner loop to calculate pairwise point distances
          a = raw_data[i]
          b = raw_data[j+i] #each simplex is a set (no order), hence [0,1] = [1,0]; so only store one
          if (i != j+i):
              dist = metric(a,b)
              if dist <= epsilon:
                  edges.append({i,j+i}) #add edge
                  weights.append([len(edges)-1,dist]) #store index and weight
  return nodes,edges,weights

# Get all the indicent edges of "node"
# where the other has lower index than "node" (to prevent replication)
def lower_nbrs(nodeSet, edgeSet, node):
    return {x for x in nodeSet if {x,node} in edgeSet and node > x}

# Calculate the complex induced by the graph, with at most (k+2) dims
# See https://blog.cc98.site/School-Courses/SP-SU-2024/TDA/lec3/#vr for further explanation
def rips(graph, k):
    nodes, edges = graph[0:2]
    VRcomplex = [{n} for n in nodes]
    for e in edges: #add 1-simplices (edges)
        VRcomplex.append(e)
    for i in range(k):
        # We are going to construct (i+2)-simplices via (i+1)-simplices
        # And since (i+1)-simplices has (i+2) points
        # We are going to use all simplices with i+2 points in this round
        for simplex in filter(lambda x: len(x)==i+2, VRcomplex): #skip 0-simplices
            # nbrs is the set of nodes, where each node is 
            # 1. incident to all nodes in the "simplex"
            # 2. has smaller index than that of all the nodes in the "simplex" (to prevent replication)
            nbrs = set.intersection(*[lower_nbrs(nodes, edges, z) for z in simplex])
            # Each node in nbrs can 
            for nbr in nbrs:
                VRcomplex.append(set.union(simplex,{nbr}))
    return VRcomplex

def drawComplex(origData, ripsComplex, axes=[-6,8,-6,6]):
  plt.clf()
  plt.axis(axes)
  plt.scatter(origData[:,0],origData[:,1]) #plotting just for clarity
  for i, txt in enumerate(origData):
      plt.annotate(i, (origData[i][0]+0.05, origData[i][1])) #add labels

  #add lines for edges
  for edge in [e for e in ripsComplex if len(e)==2]:
      #print(edge)
      pt1,pt2 = [origData[pt] for pt in [n for n in edge]]
      #plt.gca().add_line(plt.Line2D(pt1,pt2))
      line = plt.Polygon([pt1,pt2], closed=None, fill=None, edgecolor='r')
      plt.gca().add_line(line)

  #add triangles
  for triangle in [t for t in ripsComplex if len(t)==3]:
      pt1,pt2,pt3 = [origData[pt] for pt in [n for n in triangle]]
      line = plt.Polygon([pt1,pt2,pt3], closed=False, color="blue",alpha=0.3, fill=True, edgecolor=None)
      plt.gca().add_line(line)
  plt.show()
