from pyhyp import pyHyp
fileName = 'even_sphere.fmt'
#fileName = 'uneven_sphere.fmt'
#fileName = 'uneven_sphere_large.fmt'

options= {

    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile':fileName,

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 73, 
    's0':1e-4,
    'rMin':10,
 
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':1e-4,
    'pGridRatio':1.1,
    'cMax': 5,
    
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.0001,
    'volSmoothIter': 20,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 1500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('sphere.cgns')

