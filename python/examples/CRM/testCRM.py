import sys, os, time
from mdo_import_helper import import_modules, MPI, mpiPrint
import pyHyp
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65, 
    's0':1e-3,
    'gridRatio':1.20,
    'rMin':5,

    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'NMax':400,
    'ps0':1e-3,
    'pGridRatio':1.15,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 2.0,
    'epsI': 4.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.001,
    'volSmoothIter': 10,

   # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-5,
    'kspMaxIts': 500,
    'preConLag': 10,
    'kspSubspaceSize':50,
    }


hyp = pyHyp.pyHyp('3d',fileName='dpw.fmt', options=options, mirror=True)
hyp.run()
hyp.writeCGNS('dpw.cgns')
hyp.writeCGNSOrig('dpwOrig.cgns')
