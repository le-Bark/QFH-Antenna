import matplotlib.pyplot as plt
import numpy
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy

def generate_loop(diameter,height,rotation,nturns,npts):
    backfire = 1
    
    radius = diameter/2
    dir = -1
    if backfire:
        dir = 1
        height = -height

    h1_x = [0] * npts
    h1_y = [0] * npts
    h1_z = list(numpy.linspace(0,height,npts))
    h2_x = [0] * npts
    h2_y = [0] * npts
    h2_z = h1_z

    for i in range(0,npts):
        h1_x[i] = radius * math.cos(dir*i*2*math.pi*nturns/(npts-1) + rotation)
        h1_y[i] = radius * math.sin(dir*i*2*math.pi*nturns/(npts-1) + rotation)
        h2_x[i] = radius * math.cos(dir*i*2*math.pi*nturns/(npts-1) + math.pi + rotation)
        h2_y[i] = radius * math.sin(dir*i*2*math.pi*nturns/(npts-1) + math.pi + rotation)

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
loop1 = numpy.array(generate_loop(30,100,0,0.5,10))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.plot(loop1[:,0],loop1[:,1],loop1[:,2], label = "1")
plt.show()
"""