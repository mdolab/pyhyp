from pyhyp import pyHypMulti

fileName1 = 'even_sphere.fmt'
fileName2 = 'uneven_sphere.fmt'
fileName3 = 'uneven_sphere_large.fmt'

commonOptions= {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    'inputFile':fileName1,
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

# Set up specific options that will overwrite the commonOptions

options1 = {'inputFile':fileName1,
            'skip':False}
options2 = {'inputFile':fileName2,
            'skip':False}
options3 = {'inputFile':fileName3,
            'skip':False}

# Create list with different options
options = [options1,
           options2,
           options3]

# set skip list, if any
skipList = []

hyp = pyHypMulti(options=options,commonOptions=commonOptions,skipList=skipList)
hyp.combineCGNS(combinedFile='combined.cgns')
