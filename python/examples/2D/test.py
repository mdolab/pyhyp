import sys, time
sys.path.append('../../')
import pyHyp
import numpy

import petsc4py 
petsc4py.init(sys.argv)

print '2D Code is currently broken. Please check back later'
sys.exit(0)

# # Shpere with non-uniform spacing
# N = 250
# X = numpy.zeros((N,2))
# for i in xrange(N):
#     theta = 2*numpy.pi*i/(N-1) + .80*numpy.sin(2*numpy.pi*i/(N-1))
#     X[i,0] = 1*numpy.cos(theta)
#     X[i,1] = 1*numpy.sin(theta)

# # Box with 1/4 of it cut out - even spacing
# N = 255
# M = 127
# x = []
# y = []
# side = 1.0
# # Bottom
# for i in xrange(N):
#     x.append(side*i/(N-1))
#     y.append(0.0)

# # Half Right size
# for i in xrange(1,M):
#     x.append(1.0)
#     y.append(0.5*side*i/(M-1))

# # Going left
# for i in xrange(1,M):
#     x.append(1.0 - 0.5*side*i/(M-1))
#     y.append(.5)

# # Going Up
# for i in xrange(1,M):
#     x.append(.5)
#     y.append(.5 + 0.5*side*i/(M-1))

# # # Going left and up
# # for i in xrange(1,M):
# #     x.append(1.0 - 0.5*side*i/(M-1))
# #     y.append(.5 + 0.5*side*i/(M-1))

# # Going Left
# for i in xrange(1,M):
#     x.append(.5 - 0.5*side*i/(M-1))
#     y.append(1.0)

# # Going Down
# for i in xrange(1,N):
#     x.append(0.0)
#     y.append(1.0 - side*i/(N-1))
# X = numpy.array([x,y]).T

# # # -----------------------------------------
# # L domain --- looks like wing with winglet
# N = 129*1-1
# M = 65*1-1
# L = 9*1-1
# x = []
# y = []
# side = 1.0
# h = .25
# w = .025
# # Bottom
# for i in xrange(N):
#     x.append(side*i/(N-1))
#     y.append(0.0)

# # Going up
# for i in xrange(1,M):
#     x.append(side)
#     y.append(h*side*i/(M-1))

# # Going left
# for i in xrange(1,L):
#     x.append(side - w*side*i/(L-1))
#     y.append(h)

# # Going Down
# for i in xrange(1,M):
#     x.append(side - w)
#     y.append(h - (h-w)*side*i/(M-1))

# # Going left 
# for i in xrange(1,2*N):
#     x.append((side - w) - (side-w)*side*i/(2*N-1))
#     y.append(w)

# X = numpy.array([x,y]).T

# options= {'N': 97,
#           's0':5e-6,
#           'epsI':1.0,
#           'epsE':0.5,
#           'theta':3.0,
#           'volCoef':1.0,
#           'volBlend':0.01,
#           'volSmoothIter':5,
#           'gridRatio':1.15,
#           }

# options= {'N':110,
#           'Ndebug':110,
#           's0':5e-6,
#           'epsI':1.0,
#           'epsE':0.5,
#           'theta':3.0,
#           'volCoef':1.0,
#           'volBlend':0.001,
#           'volSmoothIter':0,
#           'gridRatio':1.15,
#           }

# #hyp = pyHyp.pyHyp('2d',X=X, options=options, flip=True)
# #hyp = pyHyp.pyHyp('2d',fileName='naca0012.dat', options=options,flip=True)
# hyp = pyHyp.pyHyp('2d',fileName='rae2822_marco.dat', options=options,flip=True)
# #from matplotlib.pylab import *
# #plot(X[:,0],X[:,1],'ko-')
# #show()

# timeA = time.time()
# hyp.run()
# print 'Run time:',time.time()-timeA
# hyp.writePlot3D('test.fmt')
# hyp.writeCGNS('test.cgns')
