from pyhyp import pyHyp
import numpy
options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 33,
    's0':1e-2,
    'rMin':100,

    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':.50e-7,
    'pGridRatio':1.1,
    'cMax':6.0,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
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

alpha = numpy.linspace(0,2*numpy.pi,129)
x = numpy.cos(alpha)*0.5 + 0.5
y = numpy.zeros_like(x)
for i in xrange(len(x)):
    if i < len(x)/2:
        y[i] = 0.6*(0.2969*numpy.sqrt(x[i]) - 0.1260*x[i]- 0.3516*x[i]**2 + 0.2843*x[i]**3 - 0.1036*x[i]**4)
    else:
        y[i] = -0.6*(0.2969*numpy.sqrt(x[i]) - 0.1260*x[i]- 0.3516*x[i]**2 + 0.2843*x[i]**3 - 0.1036*x[i]**4)

X = numpy.hstack([[x,y]]).T
hyp = pyHyp.pyHyp('2d',X=X, options=options, flip=True)
hyp.run()
hyp.writeCGNS('naca0012.cgns')
