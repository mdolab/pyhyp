from pyhyp import pyHyp

fileName = 'face_splay.cgns'
fileType = 'CGNS'

options= {
    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile': fileName,
    'fileType': fileType,

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65, 
    's0': 1e-6,
    'rMin': 2.5,

    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0': 1e-6,
    'pGridRatio': 1.15,
    'cMax': 5,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 0.0,
    'volCoef': 0.3,
    'volBlend': 0.001,
    'volSmoothIter': 10,

    # ---------------------------
    #   BC parameters
    # ---------------------------
    'sigmaSplay': 0.4,
    'nuSplay': 0.95,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-15,
    'kspMaxIts': 1500,
    'preConLag': 10,
    'kspSubspaceSize':50,
    'writeMetrics': False,
    }


hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('face3D.cgns')
