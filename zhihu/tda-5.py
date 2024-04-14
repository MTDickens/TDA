import numpy as np
import matplotlib.pyplot as plt

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

eps = 1.9
graph = SimplicialComplex.buildGraph(raw_data=newData, epsilon=eps)
ripsComplex = SimplicialComplex.rips(graph=graph, k=3)
SimplicialComplex.drawComplex(origData=newData, ripsComplex=ripsComplex)