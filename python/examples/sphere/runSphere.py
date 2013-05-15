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
    'N': 73, 
    's0':1e-4,
    'rMin':10,
    
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'volCoef': 0.16,
    'volBlend': 0.0005,
    'volSmoothIter': 15,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 1500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

#hyp = pyHyp.pyHyp('3d',fileName='even_sphere.fmt', options=options)
#hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere.fmt', options=options)
hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere_large.fmt', options=options)

hyp.run()
hyp.writeCGNS('sphere.cgns')
#hyp.writeCGNSOrig('sphereOrig.cgns')
