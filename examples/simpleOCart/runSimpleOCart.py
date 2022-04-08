import os
from pyhyp.utils import simpleOCart

baseDir = os.path.dirname(os.path.abspath(__file__))

# Usually this would be a nearfield volume file
# We use the corner surface file for convenience
nearFile = os.path.join(baseDir, "../corner/corner.cgns")

# Set the other inputs
dh = 0.05
hExtra = 10.0
nExtra = 49
sym = "y"
mgcycle = 3
outFile = os.path.join(baseDir, "simpleOCart.cgns")
userOptions = {"cMax": 5.0}

# Run simpleOCart
simpleOCart(nearFile, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=userOptions)
