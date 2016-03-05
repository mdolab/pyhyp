from pyhyp import pyHyp

fileName = 'plate_surf.cgns'
fileType = 'CGNS'

options= {

    # ---------------------------
    #        Input Parameters
    # ---------------------------
    'inputFile':fileName,
    'fileType': 'CGNS',
    'unattachedEdgesAreSymmetry':False,
    'outerFaceBC':'overset',
    'autoConnect':True,
    'BC':{1:{'jLow':'XYConst',
             'iLow':'XConst', 
             'iHigh':'XConst'}},
    'families':'wall',

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65, 
    's0': 1e-6,
    'marchDist':2.5,
    'splay':0.5,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0': -1,
    'pGridRatio': -1,
    'cMax': 5,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': 0.25,
    'volBlend': 0.0001,
    'volSmoothIter': 100,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('face3D.cgns')
