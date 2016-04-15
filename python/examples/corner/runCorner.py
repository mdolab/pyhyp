from pyhyp import pyHyp, pyHypMulti

'''
In this example we extrude a 90 deg corner twice using two different
setting.
This example shows how the dictionary format can be used to set up
and run multiple cases at once
'''

fileName = 'corner.cgns'
fileType = 'CGNS'

commonOptions= {

    # ---------------------------
    #        Input Parameters
    # ---------------------------
    'inputFile':fileName,
    'fileType': 'CGNS',
    'unattachedEdgesAreSymmetry':False,
    'outerFaceBC':'farfield',
    'autoConnect':True,
    'BC':{},
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

# Now set up specific options
options1 = {'outputFile':'corner1_hyp.cgns'}

options2 = {'epsE':4.0,
            'epsI':8.0,
            'outputFile':'corner2_hyp.cgns'}

# Gather options in a list
options = [options1,
           options2]

hyp = pyHypMulti(options=options,commonOptions=commonOptions)
hyp.combineCGNS()
