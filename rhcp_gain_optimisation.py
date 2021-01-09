import math
import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from generate_loop import generate_loop
import necpp

def handle_nec(result):
  if (result != 0):
    print (necpp.nec_error_message())

def antennna_height(diam,length,nturn):
    return math.sqrt((length-diam)**2 - (math.pi*diam*nturn)**2)

def add_nec_wire(nec,wire_no,p1,p2,wdiam,lmin):
    lseg = numpy.linalg.norm(p1-p2)
    nseg = math.ceil(lseg/lmin)
    handle_nec(necpp.nec_wire(nec,wire_no,nseg,p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],wdiam,1,1))

def add_nec_wire_file(f,wire_no,p1,p2,wdiam,lmin):
    lseg = numpy.linalg.norm(p1-p2)
    nseg = math.ceil(lseg/lmin)
    geometry_line = "GW {0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(wire_no,nseg,p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],wdiam)
    f.write(geometry_line)


def generate_antenna_file(freq,diam,l1,l2,nturn,wdiam):
    # l1<l2
    if l1>l2:
        print("error l1>l2")
        return -1

    l_min = l1*0.08
    
    height_diff = wdiam + 1/1000

    loop1 = numpy.array(generate_loop(diam,l1/2,0,nturn,math.ceil(l1/l_min)))
    loop2 = numpy.array(generate_loop(diam,l2/2,math.pi/2,nturn,math.ceil(l2/l_min)))

    loop1 = numpy.concatenate(([[0,0,0]],loop1,[[0,0,0]]))
    loop2 = numpy.concatenate(([[0,0,0]],loop2,[[0,0,0]]))
    #loop2 = loop2 + [0,0,height_diff] 

    f = open("out.NEC","w")
    f.write("CM test helix\n")
    f.write("CE\n")

    wireno = 0
    for i in range(len(loop1[:,0])-1):
        wireno = wireno + 1
        add_nec_wire_file(f,wireno,loop1[i,:],loop1[i+1,:],wdiam,l_min)

    wireno2 = wireno + 1

    for i in range(len(loop2[:,0])-1):
        wireno = wireno + 1
        add_nec_wire_file(f,wireno,loop2[i,:],loop2[i+1,:],wdiam,l_min)


    #add exitation card
    f.write("GE 0\n")
    f.write("EK 0\n")
    f.write("GN -1\n")
    f.write("FR  0    1    0    0    {0}       0.0\n".format(freq))
    f.write("EX 5 {0} {1} 0 1.0 0.0\n".format(1,1))
    f.write("EX 5 {0} {1} 0 -1.0 0.0\n".format(wireno2,1))
    f.write("RP  0   35   35 1001       0.0       0.0      10.0      10.0\n")
    f.write("EN\n")

    f.close()



def antenna_freq_annalysis(freq,diam,l1,l2,nturn,wdiam):

    # l1<l2
    if l1>l2:
        print("error l1>l2")
        return 0

    l_min = l1*0.08
    
    height_diff = wdiam + 1/1000

    loop1 = numpy.array(generate_loop(diam,l1/2,0,nturn,math.ceil(l1/l_min)))
    loop2 = numpy.array(generate_loop(diam,l2/2,math.pi/2,nturn,math.ceil(l2/l_min)))

    loop1 = numpy.concatenate(([[0,0,0]],loop1,[[0,0,0]]))
    loop2 = numpy.concatenate(([[0,0,0]],loop2,[[0,0,0]]))
    loop2 = loop2 + [0,0,height_diff] 

    nec = necpp.nec_create()

    wireno = 0
    for i in range(len(loop1[:,0])-1):
        wireno = wireno + 1
        add_nec_wire(nec,wireno,loop1[i,:],loop1[i+1,:],wdiam,l_min)

    wireno2 = wireno + 1

    for i in range(len(loop2[:,0])-1):
        wireno = wireno + 1
        add_nec_wire(nec,wireno,loop2[i,:],loop2[i+1,:],wdiam,l_min)

    handle_nec(necpp.nec_geometry_complete(nec,0))

    #ground card, no ground
    handle_nec(necpp.nec_gn_card(nec,-1,0,0,0,0,0,0,0))
    #extended wire kernel
    handle_nec(necpp.nec_ek_card(nec,0))

    handle_nec(necpp.nec_ex_card(nec,5,1,1,0,1,0,0,0,0,0))
    handle_nec(necpp.nec_ex_card(nec,5,wireno2,1,0,-1,0,0,0,0,0))

    handle_nec(necpp.nec_fr_card(nec,0,1,freq,0))
    handle_nec(necpp.nec_rp_card(nec,0,35,35,1,0,0,1,0,0,10,10,0,0))
    
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    plt.plot(loop1[:,0],loop1[:,1],loop1[:,2], label = "1")
    plt.plot(loop2[:,0],loop2[:,1],loop2[:,2], label = "2")
    plt.show()
    """
    return nec

def max_rhcp_gain_of_antenna(freq,diam,l1,l2,nturn,wdiam):
    # l1<l2
    if l1>l2:
        print("error l1>l2")
        return 0
    nec = antenna_freq_annalysis(freq,diam,l1,l2,nturn,wdiam)
    gain = necpp.nec_gain_rhcp_max(nec,0)
    impedance = complex(necpp.nec_impedance_real(nec,0),necpp.nec_impedance_imag(nec,0))
    necpp.nec_delete(nec)
    return gain

def impedance_of_antenna(freq,diam,l1,l2,nturn,wdiam):
    # l1<l2
    if l1>l2:
        print("error l1>l2")
        return -1
    nec = antenna_freq_annalysis(freq,diam,l1,l2,nturn,wdiam)
    gain = necpp.nec_gain_rhcp_max(nec,0)
    impedance = complex(necpp.nec_impedance_real(nec,0),necpp.nec_impedance_imag(nec,0))
    necpp.nec_delete(nec)
    return impedance


diameter = 0.03
def f(x):
    print(x)
    res = max_rhcp_gain_of_antenna(1575,diameter,x[0],x[1],0.5,0.001)
    print(res)
    return -1*res

def constraint1(x):
    return x[1]-x[0]
def constraint2(x):
    return impedance_of_antenna(1575,diameter,x[0],x[1],0.5,0.001).real-50
def constraint3(x):
    return impedance_of_antenna(1575,diameter,x[0],x[1],0.5,0.001).imag

x0 = [0.190,0.200]

from scipy.optimize import minimize

bnds = ((0.190,0.23),(0.190,0.23))
con1 = {'type': 'ineq', 'fun': constraint1}
con2 = {'type': 'eq', 'fun': constraint2}
con3 = {'type': 'eq', 'fun': constraint3}
cons = [con1]
#print(max_rhcp_gain_of_antenna(1575,diameter,0.190,0.191,0.5,1/1000))

sol = minimize(f,x0,method='SLSQP',bounds=bnds,constraints=cons,options={'disp': 'true'})
print(sol)
print("length 1 : {0}".format(sol.x[0]))
print("height 1 : {0}".format(antennna_height(diameter,sol.x[0]/2,0.5)))
print("length 2 : {0}".format(sol.x[1]))
print("height 2 : {0}".format(antennna_height(diameter,sol.x[1]/2,0.5)))

generate_antenna_file(1575,diameter,sol.x[0],sol.x[1],0.5,1/1000)
