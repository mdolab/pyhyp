from pyhyp import pyHyp
from pygeo import pyBlock, DVGeometry
from pyspline import pySpline
import numpy
# Define params here:
ffd_file = 'mdo_tutorial_ffd.fmt'

options= {
   # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65,#145,#73, 
    's0':5e-6, #5e-6
    'rMin':23.2,

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
    'kspRelTol': 1e-5,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

def winglet(val, geo):

    # Move the last point on ref axis:

    C = geo.extractCoef(0)
    s = geo.refAxis.curves[0].s

    C[-1,0] += val[0]
    C[-1,1] += val[1]
    C[-1,2] += val[2]

    # Also need to get "dihedreal" angle for this section
    theta = numpy.arctan(val[1]/(C[-1,2] - C[-2,2]))
    geo.rot_x[0].coef[-1] = -theta*180/numpy.pi

    geo.restoreCoef(C, 0)

    # Set the chord as well
    geo.scale[0].coef[-1] = val[3]

    geo.restoreCoef(C, 0)

    return

# -> Define the volume ID's that coorespond to the wing and tail:
wing_vols = [0]
FFD = pyBlock.pyBlock('plot3d', file_name=ffd_file, 
                      file_type='ascii', order='f', FFD=True)

coef = FFD.vols[0].coef.copy()
nSpan = coef.shape[2]
ref = numpy.zeros((nSpan,3))

for k in xrange(nSpan):
    max_x = numpy.max(coef[:,:,k,0])
    min_x = numpy.min(coef[:,:,k,0])

    ref[k,0] = min_x + 0.25*(max_x-min_x)
    ref[k,1] = numpy.average(coef[:,:,k,1])
    ref[k,2] = numpy.average(coef[:,:,k,2])

c0 = pySpline.curve(X=ref,k=2)

# Create the hyp object
#hyp = pyHyp.pyHyp('3d',fileName='717_large.fmt', options=options, zMirror=True)
hyp = pyHyp.pyHyp('3d',fileName='717_small.fmt', options=options, zMirror=True)

coords = hyp.X.copy()
nx = len(hyp.X)
ind = []
for i in xrange(len(hyp.X)):
    if coords[i][2] < 0:
        coords[i][2] = -coords[i][2]
        ind.append(i)

DVGeo = DVGeometry.DVGeometry([coords], [c0], vol_list=[wing_vols],
                              axis = ['x'],FFD=FFD,rot_type=5,
                              names=['warp_coords'])


DVGeo.addGeoDVGlobal('winglet',[0,0,0,1],-5,5,winglet)
DVGeo.setValues('winglet',[1.5,2.5,-2.0,.60],scaled=False) # Large, good

newCoords = DVGeo.update('warp_coords')
for i in xrange(len(ind)):
    newCoords[ind[i]][2] = -newCoords[ind[i]][2]

# Set the 'new' coordinates in the hypObject. If you want the
# undeflected wing, comment out the line below.
hyp.X = newCoords
hyp.run()
hyp.writeCGNS('717.cgns')
#hyp.writeCGNSOrig('717Orig.cgns')
