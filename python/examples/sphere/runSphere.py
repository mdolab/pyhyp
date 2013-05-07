import sys, os, time
from mdo_import_helper import import_modules, MPI, mpiPrint
sys.path.append('../../')
import pyHyp
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 32, 
    's0':1e-2,
    'rMin':25,
    'gridRatio':1.2,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 0.5,
    'volCoef': 0.16,
    'volBlend': 0.001,
    'volSmoothIter': 5,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-5,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

#hyp = pyHyp.pyHyp('3d',fileName='even_sphere.fmt', options=options)
#hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere.fmt', options=options)
hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere_large.fmt', options=options)

hyp.run()
hyp.writeCGNS('sphere.cgns')
#hyp.writeCGNSOrig('sphereOrig.cgns')
