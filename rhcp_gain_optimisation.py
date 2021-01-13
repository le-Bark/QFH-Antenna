import math
import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from generate_loop import generate_loop
import necpp
from scipy.optimize import minimize

def handle_nec(result):
  if (result != 0):
    print (necpp.nec_error_message())
  return result

def antennna_height(diam,length,nturn):
    return math.sqrt((length-diam)**2 - (math.pi*diam*nturn)**2)

def add_nec_wire(nec,wire_no,p1,p2,wdiam,lmin):
    lseg = numpy.linalg.norm(p1-p2)
    nseg = math.ceil(lseg/lmin)
    return handle_nec(necpp.nec_wire(nec,wire_no,nseg,p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],wdiam,1,1))

def add_nec_wire_file(f,wire_no,p1,p2,wdiam,lmin):
    lseg = numpy.linalg.norm(p1-p2)
    nseg = math.ceil(lseg/lmin)
    geometry_line = "GW {0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(wire_no,nseg,p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],wdiam)
    f.write(geometry_line)


def generate_antenna_file(freq,diam1,diam2,h1,h2,nturn,wdiam):

    l_min = h1*0.08
    
    height_diff = wdiam + 1/1000

    loop1 = numpy.array(generate_loop(diam1,h1,0,nturn,math.ceil(h1/l_min)))
    loop2 = numpy.array(generate_loop(diam2,h2,math.pi/2,nturn,math.ceil(h2/l_min)))

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



def antenna_freq_annalysis(freq,diam1,diam2,h1,h2,nturn,wdiam):

    l_min = h1*0.08
    
    height_diff = wdiam + 1/1000

    loop1 = numpy.array(generate_loop(diam1,h1,0,nturn,math.ceil(h1/l_min)))
    loop2 = numpy.array(generate_loop(diam2,h2,math.pi/2,nturn,math.ceil(h2/l_min)))

    loop1 = numpy.concatenate(([[0,0,0]],loop1,[[0,0,0]]))
    loop2 = numpy.concatenate(([[0,0,0]],loop2,[[0,0,0]]))
    loop2 = loop2 + [0,0,height_diff] 
    #l1>l2

    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    plt.plot(loop1[:,0],loop1[:,1],loop1[:,2], label = "1")
    plt.plot(loop2[:,0],loop2[:,1],loop2[:,2], label = "2")
    plt.show()
    """

    nec = necpp.nec_create()
    necerror = 0
    wireno = 0
    for i in range(len(loop1[:,0])-1):
        wireno = wireno + 1
        necerror = necerror + add_nec_wire(nec,wireno,loop1[i,:],loop1[i+1,:],wdiam,l_min)

    wireno2 = wireno + 1

    for i in range(len(loop2[:,0])-1):
        wireno = wireno + 1
        necerror = necerror + add_nec_wire(nec,wireno,loop2[i,:],loop2[i+1,:],wdiam,l_min)



    if necerror:
        return 0

    handle_nec(necpp.nec_geometry_complete(nec,0))

    #ground card, no ground
    handle_nec(necpp.nec_gn_card(nec,-1,0,0,0,0,0,0,0))
    #extended wire kernel
    handle_nec(necpp.nec_ek_card(nec,0))

    handle_nec(necpp.nec_ex_card(nec,5,1,1,0,1,0,0,0,0,0))
    handle_nec(necpp.nec_ex_card(nec,5,wireno2,1,0,-1,0,0,0,0,0))

    handle_nec(necpp.nec_fr_card(nec,0,1,freq,0))
    handle_nec(necpp.nec_rp_card(nec,0,35,35,1,0,0,1,0,0,10,10,0,0))
    
    
    return nec

def max_rhcp_gain_of_antenna(freq,diam1,diam2,l1,l2,nturn,wdiam):

    nec = antenna_freq_annalysis(freq,diam1,diam2,l1,l2,nturn,wdiam)
    if nec == 0:
        return -10
    gain = necpp.nec_gain_rhcp_max(nec,0)
    impedance = complex(necpp.nec_impedance_real(nec,0),necpp.nec_impedance_imag(nec,0))
    necpp.nec_delete(nec)
    return gain

def return_loss_of_antenna(freq,diam1,diam2,l1,l2,nturn,wdiam):
    
    nec = antenna_freq_annalysis(freq,diam1,diam2,l1,l2,nturn,wdiam)
    if nec == 0:
        return -10
    gain = necpp.nec_gain_rhcp_max(nec,0)
    impedance = complex(necpp.nec_impedance_real(nec,0),necpp.nec_impedance_imag(nec,0))
    return_loss = 20*math.log10(abs((impedance-50)/(impedance+50)))
    necpp.nec_delete(nec)
    return return_loss

freqmhz = 1575
wire_diameter = 0.001
wh_ratio = 0.25
wavelength = 3e8 / (freqmhz * 1e6)
nturns = 0.5

# diameter height
x0 = [wavelength/4,wavelength/4]
bnds = ((wavelength/10,wavelength/2),(wavelength/10,wavelength/2))

def f1(x):
    print(x)
    res = return_loss_of_antenna(freqmhz,x[0],x[0],x[1],x[1],nturns,wire_diameter)
    print(res)
    return res

#w/h constraint
def constraint1(x):
    return x[0]/x[1]-wh_ratio
def constraint2(x):
    return math.sqrt((x[0]*math.pi*nturns)**2 + (x[1])**2)+x[0]-wavelength/2

con1 = {'type': 'eq', 'fun': constraint1}
con2 = {'type': 'ineq', 'fun': constraint2}
cons = [con1,con2]

sol = minimize(f1,x0,method='SLSQP',bounds=bnds,constraints=cons,options={'disp': 'true'})
print(sol)


diam1 = sol.x[0]
height1 = sol.x[1]

def f2(x):
    print(x)
    res = max_rhcp_gain_of_antenna(freqmhz,x[0],x[1],x[2],x[3],nturns,wire_diameter)
    print(res)
    return -res

bnds = ((diam1,1.2*diam1),(0.8*diam1,diam1),(height1,1.2*height1),(0.8*height1,height1))
x0 = [1.01*diam1,0.99*diam1,1.01*height1,0.99*height1]

sol = minimize(f2,x0,method='SLSQP',bounds=bnds,options={'disp': 'true'})

print(sol)
print("calculated wavelength : {0}".format((math.sqrt((diam1*math.pi*nturns)**2 + (height1)**2)+diam1)*2))
print("diam1 : {0}".format(diam1))
print("height1 : {0}".format(height1))
generate_antenna_file(freqmhz,sol.x[0],sol.x[1],sol.x[2],sol.x[3],nturns,wire_diameter)