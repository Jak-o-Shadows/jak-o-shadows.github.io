# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:10:12 2015

@author: Jak
"""

import math
import numpy


#define coord system origin as the centre of the bottom plate

#Find base plate attachment locations
bAngles = [15, 105, 135, 225, 255, 345]
bAngles = [math.radians(x) for x in bAngles]
bR = 50
bPos = [[bR*math.cos(theta), bR*math.sin(theta), 0] for theta in bAngles]
bPos = numpy.array(bPos)

#Platform attachment locations
pAngles = [45, 75, 165, 195, 285, 315]
pAngles = [math.radians(x) for x in pAngles]
pR = 50
pPos = [[pR*math.cos(theta), pR*math.sin(theta), 0] for theta in pAngles]
pPos = numpy.array(pPos)


#Input lengths of sensors
L = numpy.array([122.759, 122.759, 122.759, 122.759, 122.759, 122.759]).transpose()


#newton-raphson
tol_f = 1e-3;
tol_a = 1e-3;
#iteration limits
maxIters = 1e3
iterNum = 0

#initial guess position
# a = [x, y, z, phi, theta, psi] - angles in degrees initially
a = [20, 20, 100, 10, 10, 10]
a[3:] = [math.radians(x) for x in a[3:]] #convert to radians
a = numpy.array(a).transpose()
while iterNum < maxIters:
    iterNum += 1
    
    x = a[0]
    y = a[1]
    z = a[2]
    phi = a[3]
    th = a[4]
    psi = a[5]
    #Must translate platform coordinates into base coordinate system
    #Calculate rotation matrix elements
    cphi = math.cos(phi)
    sphi = math.sin(phi)
    cth = math.cos(th)
    sth = math.sin(th)
    cpsi = math.cos(psi)
    spsi = math.sin(psi)   
    #Hence calculate rotation matrix
    #Note that it is a 3-2-1 rotation matrix
    Rzyx = numpy.array([[cpsi*cth, cpsi*sth*sphi - spsi*cphi, cpsi*sth*cphi + spsi*sphi] \
                        ,[spsi*cth, spsi*sth*sphi + cpsi*cphi, spsi*sth*cphi - cpsi*sphi] \
                        , [-sth, cth*sphi, cth*cphi]])
                        
    #Hence platform sensor points with respect to the base coordinate system
    xbar = a[0:3] - bPos

    #Hence orientation of platform wrt base

    uvw = numpy.zeros(pPos.shape)
    for i in xrange(6):
        uvw[i, :] = numpy.dot(Rzyx, pPos[i, :])

    
    #Hence find value of objective function
    #The calculated lengths minus the actual length
    f = -1 * (numpy.sum(numpy.square(xbar + uvw), 1) - numpy.square(L))
    sumF = numpy.sum(numpy.abs(f))
    if sumF < tol_f:
        #success!
        print "Converged! Output is in 'a' variable"
        break
    
    #As using the newton-raphson matrix, need the jacobian (/hessian?) matrix
    #Using paper linked above:
    dfda = numpy.zeros((6, 6))
    dfda[:, 0:3] = 2*(xbar + uvw)
    for i in xrange(6):
        #Numpy * is elementwise multiplication!!
        #Indicing starts at 0!
        #dfda4 is swapped with dfda6 for magic reasons!  
        dfda[i, 5] = 2*(-xbar[i,0]*uvw[i,1] + xbar[i,1]*uvw[i,0]) #dfda4
        dfda[i, 4] = 2*((-xbar[i,0]*cpsi + xbar[i,1]*spsi)*uvw[i,2] \
                        - (pPos[i,0]*cth + pPos[i,1]*sth*sphi)*xbar[i,2]) #dfda5
        dfda[i, 3] = 2*pPos[i, 1]*(numpy.dot(xbar[i,:],Rzyx[:,2])) #dfda

    #Hence solve system for delta_{a} - The change in lengths
    delta_a = numpy.linalg.solve(dfda, f)

    if abs(numpy.sum(delta_a)) < tol_a:
        print "Small change in lengths -- converged?"
        break
    a = a + delta_a

for i in xrange(3,6):
    a[i] = math.degrees(a[i])
print a
print iterNum

