"""
This script uses the NACA 0012 airfoil equation to generate the 
datapoints and a 2D RANS mesh. This mesh has a blunt trailing edge.
"""


from pyhyp import pyHyp
import numpy
import matplotlib.pyplot as plt
options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 129,
    's0':0.000001,
    'rMin':100,

    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':1e-8,
    'pGridRatio':1.1,
    'cMax':6.0,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 0, 
    'epsI': 2.0,
    'theta': 4.0,
    'volCoef': .16,
    'volBlend': 0.005,
    'volSmoothIter': 25,
    
    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-15,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

alpha = numpy.linspace(0, 2*numpy.pi, 273)
x = numpy.cos(alpha)*0.5 + 0.5
y = numpy.zeros_like(x)

for i in xrange(len(x)):
    if i < len(x)/2:
        y[i] = 0.6*(0.2969*numpy.sqrt(x[i]) - 0.1260*x[i] - 0.3516*x[i]**2 + 0.2843*x[i]**3 - 0.1015*x[i]**4)
    else:
        y[i] = -0.6*(0.2969*numpy.sqrt(x[i]) - 0.1260*x[i] - 0.3516*x[i]**2 + 0.2843*x[i]**3 - 0.1015*x[i]**4)

# Since the TE is open we need to close it. Close it multiple linear segments.
delta_y = numpy.linspace(y[-1], y[0], 16, endpoint=True)
delta_y = delta_y[1:]

x = numpy.append(x, numpy.ones_like(delta_y))
y = numpy.append(y, delta_y)

# Combine into one array
X = numpy.hstack([[x,y]]).T

# Write out the points in case we want to plot it in tecplot
numpy.savetxt("naca0012_rans.dat",X)

# Generate the mesh
hyp = pyHyp.pyHyp('2d',X=X, options=options, flip=True)
hyp.run()
hyp.writeCGNS('naca0012_rans.cgns')
