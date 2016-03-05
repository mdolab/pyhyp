from pyhyp import pyHyp
#fileName = 'even_sphere.fmt'
#fileName = 'uneven_sphere.fmt'
fileName = 'uneven_sphere_large.fmt'

options= {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    'inputFile':fileName,
    'unattachedEdgesAreSymmetry':True,
    'outerFaceBC':'farField',
    'autoConnect':True,
    'BC':{},
    'families':'wall',

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 73, 
    's0':1e-4,
    'marchDist':10,
 
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':-1,
    'pGridRatio':-1,
    'cMax': 5,
    
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .25,
    'volBlend': 0.0001,
    'volSmoothIter': 100,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('sphere.cgns')

