import matplotlib.pyplot as plt
import numpy
from mpl_toolkits.mplot3d import Axes3D
import math


def generate_loop(diameter,length,rotation,nturns,npts):
    
    radius = diameter/2
    height = 1* math.sqrt(math.pow(length-diameter,2) - math.pow(2*math.pi*radius*nturns,2))
    
    wh = diameter/height

    #print(wh)

    h1_x = [0] * npts
    h1_y = [0] * npts
    h1_z = list(numpy.linspace(0,height,npts))
    h2_x = [0] * npts
    h2_y = [0] * npts
    h2_z = h1_z

    for i in range(0,npts):
        h1_x[i] = radius * math.cos(-i*2*math.pi*nturns/(npts-1) + rotation)
        h1_y[i] = radius * math.sin(-i*2*math.pi*nturns/(npts-1) + rotation)
        h2_x[i] = radius * math.cos(-i*2*math.pi*nturns/(npts-1) + math.pi + rotation)
        h2_y[i] = radius * math.sin(-i*2*math.pi*nturns/(npts-1) + math.pi + rotation)

    h1_x = h1_x + [0] + list(reversed(h2_x))
    h1_y = h1_y + [0] + list(reversed(h2_y))
    h1_z = h1_z + [height] + list(reversed(h2_z))

    h = [0] * len(h1_x)

    """
    h1_x.insert(0,0)
    h1_y.insert(0,0)
    h1_z.insert(0,0)
    """
    
    for i in range(len(h1_x)):
        h[i] = [h1_x[i],h1_y[i],h1_z[i]]

    return h

"""
t = generate_loop(30,100,0,0.5,10)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.plot(t[0],t[1],t[2], label = "1")
plt.show()
"""