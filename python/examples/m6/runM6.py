from pyhyp import pyHyp
fileName = 'm6_small.fmt'
#fileName = 'm6_large.fmt'
#fileName = 'm6_tiny.fmt'

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
    'N': 81,
    's0':1.5e-5,
    'marchDist':30,
    'nConstantStart':1,
    
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':-1,
    'pGridRatio':-1,
    'cMax':5.0,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .25,
    'volBlend': 0.0005,
    'volSmoothIter': 100,
    'kspreltol':1e-4,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('m6.cgns')


