from pyhyp import pyHyp
from pygeo import pyBlock, DVGeometry
from pyspline import pySpline
import numpy

# Define params here:
ffd_file = 'mdo_tutorial_ffd.fmt'
fileName = '717_small.fmt'
#fileName = '717_large.fmt'

options= {

    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile':fileName,

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65,
    's0':5e-6,
    'rMin':23.2,
    'mirror':'z',

    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':5e-6, #5e-6,
    'pGridRatio':1.1,
    'cMax': 5,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .20,
    'volBlend': 0.0005,
    'volSmoothIter': 20,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

def winglet(val, geo):

    # Move the last point on ref axis:

    C = geo.extractCoef('wing')
    s = geo.extractS('wing')

    C[-1,0] += val[0]
    C[-1,1] += val[1]
    C[-1,2] += val[2]

    # Also need to get "dihedreal" angle for this section
    theta = numpy.arctan(val[1]/(C[-1,2] - C[-2,2]))
    geo.rot_x['wing'].coef[-1] = -theta*180/numpy.pi

    geo.restoreCoef(C, 'wing')

    # Set the chord as well
    geo.scale['wing'].coef[-1] = val[3]

# Create the hyp object
hyp = pyHyp(options=options)

coords = hyp.getSurfaceCoordinates()
nx = len(coords)
ind = []
for i in xrange(nx):
    if coords[i][2] < 0:
        coords[i][2] = -coords[i][2]
        ind.append(i)

DVGeo = DVGeometry(ffd_file)
coef = DVGeo.FFD.vols[0].coef.copy()

# First determine the reference chord lengths:                                                                                                                                        
nSpan = coef.shape[2]
ref = numpy.zeros((nSpan,3))

for k in xrange(nSpan):
    max_x = numpy.max(coef[:,:,k,0])
    min_x = numpy.min(coef[:,:,k,0])

    ref[k,0] = min_x + 0.25*(max_x-min_x)
    ref[k,1] = numpy.average(coef[:,:,k,1])
    ref[k,2] = numpy.average(coef[:,:,k,2])

c0 = pySpline.Curve(X=ref,k=2)
DVGeo.addRefAxis('wing', c0)

DVGeo.addGeoDVGlobal('winglet',[0,0,0,1], winglet, lower=-5, upper=5)
DVGeo.setDesignVars({'winglet': [1.5,2.5,-2.0,.60]})

DVGeo.addPointSet(coords, 'coords')
newCoords = DVGeo.update('coords')

for i in xrange(len(ind)):
    newCoords[ind[i]][2] = -newCoords[ind[i]][2]

# Set the 'new' coordinates in the hypObject. If you want the
# undeflected wing, comment out the line below.
hyp.setSurfaceCoordinates(newCoords)
hyp.run()
hyp.writeCGNS('717.cgns')
