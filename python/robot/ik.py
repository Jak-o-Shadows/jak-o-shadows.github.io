# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:33:56 2013

@author: Jak
"""

############################################################################################
###                                                                                      ###
### PyGame with three lanes of cars coming at you and you have to dodge them! Have fun!  ###
###                                                                                      ###
### Author: SUHAS SG                                                                     ###
### jargnar@gmail.com                                                                    ###
###                                                                                      ###
### suhased.wordpress.com                                                                ###
### twitter: @jargnar                                                                    ### 
### facebook: facebook.com/jargnar                                                       ###
###                                                                                      ###
###                                                                                      ###
### Disclaimer: The kelvinized font and the small car pictures                           ###
### is not mine and I found it on the internet.                                          ###
### All the other images are mine.                                                       ###
###                                                                                      ###
### Do Enjoy the game!                                                                   ###
### You need to have Python and PyGame installed to run it.                              ###
###                                                                                      ###
### Run it by typing "python dodgecars.py" in the terminal                               ###
###                                                                                      ###
###                                                                                      ###
############################################################################################
###QUITE QUITE HEAVILY MODIFIED by Jak_o_Shadows
import pygame,sys
import math
import numpy

#global init
pygame.init()
size = 1280,1024
textcolor = 233,230,20
speed = -3
up = right = True 
down = left = False
screen = pygame.display.set_mode(size)
bg = pygame.image.load("res/bg.png")
bgrect = bg.get_rect()
pygame.key.set_repeat(65,65)

timestep = 10/100.0


class Arm:
    'A whole robot arm!'
    def __init__(self, method, screen):
        #display stuff
        size = screen.get_size()
        self.imgR = pygame.image.load("res/car3.png")
        self.img = self.imgR
        self.rect = self.imgR.get_rect()
        self.s = pygame.Surface(size, flags=0)
        self.screen = screen
        self.rect = self.s.get_rect()
        self.rect.topleft = (0, 0)
        self.center = numpy.zeros((3, 1))
        self.center[0,0] = size[0]/2
        self.center[1,0] = size[1]/2
        self.colours = {"line": (0, 255, 0, 0), "blob" : (255, 0, 0, 0), "target": (0, 0, 255,0), \
                        "shell" : (125,125,125,0)}
        
        #arm stuff
        
        self.method = "ccd"
        self.methodList = {"ccd": self.ccd, "ikT": self.ikT, "ikPinv": self.ikPinv}
        
        pos = numpy.zeros((3, 1))
        self.sections= []
        inf = float("inf")
        self.sections.append({"length":200, "start":pos, "theta": math.pi/4, "up":math.pi, "low":0})
        self.sections.append({"length":200, "start":numpy.copy(pos), "theta" : math.pi/2, "up":math.pi/2, "low":-math.pi/2})
        self.sections.append({"length":80, "start":numpy.copy(pos), "theta" : math.pi/2, "up":math.pi/2, "low":-math.pi/2})
        
        
        #graphical representation of target
        self.target = (0, 0, 0)

        #draw it
        self.drawArm()
        
    def cossin(self, h, theta):
        """simple x, y, from a radius and angle"""
        """Could I use imaginary numbers?"""
        dirMat = numpy.array([[math.cos(theta)], [math.sin(theta)], [0]])
        return h*dirMat #not matrix multiplication, just elementwise
        
    def length(self):
        """Returns the total length of the arm"""
        length = sum([x["length"] for x in self.sections])
        return length
        
    def alpha(self, jac, goal):
        error = 1
        numa = 1
        
    def inverseK(self,goal):
        self.methodList[self.method](goal)
        
    def ikPinv(self, pos, iterLimit=500):
        self.ik(pos, iterLimit, method="pinv")
        
    def ikT(self, pos, iterLimit=500):
        self.ik(pos, iterLimit, method="t")
            
    def ik(self, pos, iterLimit=500, method="t"):
        """Uses the jacobian tranverse/inverse ("t", "pinv") method for ik"""
        goal = numpy.array(pos + (0,))
        goal.shape = (3, 1)
        goal = self.flipY(goal)
        iters = 0
        errorLimit = 5**2
        coords = self.sectionsGlobalPos()
        error = self.distSquare(goal, coords[-1][1])
        while error > errorLimit:
            #find the jacobian
            jac, coords = self.jacobian()
        
            #ik time!
            dif = goal - coords[-1][1]
            if method == "t":
                a = 1e-5
                dTheta = a*numpy.dot(jac.T,dif)
            elif method=="pinv":
                dTheta = numpy.dot(numpy.linalg.pinv(jac), dif)
            else:
                print "no such method for jacobian ik"
            for index in xrange(len(self.sections)):
                self.sections[index]["theta"] += dTheta[index, 0]
                up = self.sections[index]["up"]
                low = self.sections[index]["low"]
                self.sections[index]["theta"] = self.twoSideLim(self.sections[index]["theta"], up, low)
                #limit to +- 90 joint angles
                #self.sections[i]["theta"] = self.basicLimit(self.sections[i]["theta"])
                
            #update the error
            #possibly may run one iteration over. Who cares.
            error = self.distSquare(coords[-1][1], goal)
            iters += 1
            #we want to bail out if we're not reaching it very quickly
            if iters > iterLimit:
                error = 0
                
        self.drawArm()

        
    def jacobian(self):
        """Calculates the Jacobian for the jacobian tranverse/inverse ik methods.
        Uses a first order taylor series (linear) for the approximation of each
        partial derivative
        """
        v = numpy.array([0,0,1]) #as all are rotating about the same axis (z?)
        coords = self.sectionsGlobalPos()

        difs = numpy.zeros((len(self.sections), len(self.sections)),dtype=numpy.ndarray)
        jac = numpy.zeros((3, len(self.sections)),dtype=numpy.ndarray)
        for i in xrange(len(coords)):
            for j in xrange(len(coords)):
                difs[i, j] = coords[j][1] - coords[i][0]
              
        for i in xrange(difs.shape[1]):
            q = difs[0,i]
            jac[:,i] = numpy.cross(v.T, q.T)

        
        return jac, coords
        
        
    def ccd(self, pos, iterLimit = 200, errorLimit=20, thetaStep=None):
        """Uses the Cyclic Coordinate Descent method for ik
        
        Translated from https://sites.google.com/site/auraliusproject/ccd-algorithm        
        """
        goal = self.tupToArr(pos)
        goal = self.flipY(goal)
        coords = self.sectionsGlobalPos()
        errorLimit = errorLimit**2
        error = self.distSquare(coords[-1][0], goal)
        count = 0
        while error > errorLimit:
            for linkNo in range(len(self.sections))[::-1]: #reverse order for ccd
                pe = coords[-1][1] #end effector end location
                pc = coords[linkNo][0]  
                
                a = (pe-pc)/self.dist(pe, pc)
                
                b = (goal-pc)/self.dist(goal,pc)
                
                cosTeta = numpy.dot(a.T[0], b.T[0])
                #because of probably numerical errors, limit to +- 1
                cosTeta = self.basicLimit(cosTeta, 1)
                teta = math.acos(cosTeta)
                if thetaStep:
                    teta = self.basicLimit(teta, thetaStep)
                
                dirVector = numpy.cross(a.T[0], b.T[0])
                
                #need to check what plane the link is about, because 3D
                #planeNorm = numpy.cross(                
                
                #check if ccw or cw, as cos is an even function, and doesn't give ccw or cw
                if dirVector[2] <0:
                    teta *= -1
                
                #update the current links angle
                self.sections[linkNo]["theta"] += teta
                self.sections[linkNo]["theta"] = self.resetTheta(self.sections[linkNo]["theta"])
                up = self.sections[linkNo]["up"]
                low = self.sections[linkNo]["low"]
                self.sections[linkNo]["theta"] = self.twoSideLim(self.sections[linkNo]["theta"], up, low)
                #update all the links position
                coords = self.sectionsGlobalPos()
                
            error = self.distSquare(coords[-1][0], goal)
            count += 1
            
            #we want to bail out if we're not reaching it very quickly
            if count > iterLimit:
                error = 0
        self.drawArm()
        
    def basicLimit(self, theta, lim=math.pi):
        """Limit the value of theta to between +- lim"""
        if theta>lim:
            return lim
        elif theta<-1*lim:
            return -1*lim
        else:
            return theta
            
            
    def twoSideLim(self, x, up, low):
        """Clamps the value of x to between up and low inclusive"""
        if x < up:
            if x >= low:
                return x
            else:
                return low
        else:
            return up
            
    def resetTheta(self, theta):
        """sin, cos are 2pi period, stop the arm angles from getting silly
        by limiting to +- pi (?)"""
        num = theta/(2*math.pi)
        if num > 0:
            num = math.floor(num)
        else:
            num = math.ceil(num)
        return theta - num*2*math.pi
            


    #coordinate stuff
    def trans(self, pos, origin):
        """Translates pos local to origin to parent coords"""
        return pos + origin
        
    def transR(self, theta, thetaO):
        """Translates local rotation angle to parent rotation system"""
        return theta + thetaO
        
    def rotMat(self, pos, theta):
        """Rotates pos theta around implicit local origin using rotation matrix"""
        r1 = [math.cos(theta), -math.sin(theta), 0]
        r2 = [math.sin(theta), math.cos(theta), 0]
        r3 = [0, 0, 1]
        base = [r1, r2, r3]
        mat = numpy.array(base)
        return numpy.dot(mat, pos)
        
    def dist(self, start, end):
        """Calculates the distance between the start and end points"""
        distance = math.sqrt(sum((end-start)**2))
        return distance
        
    def distSquare(self, start, end):
        """Calculates the squared distance between the start and end points"""
        return sum((end-start)**2)
        
    def tupToArr(self, pos):
        """Converts a tuple of x, y, z coordinates to the numpy array used everywhere"""
        if len(pos) == 2:
            goal = numpy.array(pos + (0,))
        else:
            goal = numpy.array(pos)
        goal.shape = (3, 1)
        return goal
        
    def sectionsGlobalPos(self):
        """Returns self.sections, but using global coordinates Also in far simpler form
        of [start, end,  theta] for each arm section"""
        coords = []
        for i in xrange(len(self.sections)): 
            start = self.sections[i]["start"]
            theta = self.sections[i]["theta"]
            l = self.sections[i]["length"]
            end = start + self.cossin(l, theta)
            if i==0:
                tstart = self.trans(start, self.center)
                tend = self.trans(end, self.center)
                ttheta = self.transR(theta, 0)
            else:
                tstart = self.trans(start, coords[-1][1])
                tend = self.trans(self.rotMat(end, coords[-1][2]), coords[-1][1])
                ttheta = self.transR(theta, coords[-1][2])
            coords.append([tstart,tend, ttheta])
        return coords
        
    
    #display stuff
    def drawArm(self):
        """Draws the arm on screen"""
        #not doing this well at all:
        self.s = pygame.Surface(size, flags=0)
        #draw the max region/accessible area of the arm
        pygame.draw.circle(self.s, self.colours["shell"], self.ints(self.center), self.length(), 0)
        #having a blank surface, we draw on it
        coords = self.sectionsGlobalPos()
        for start, end, theta in coords:
            self.line(start, end)
            self.blob(end)
        pygame.draw.circle(self.s, self.colours["target"], self.target[:2], 5, 0)
        pygame.draw.circle(self.s, (255, 255, 255, 0), self.ints(self.center), 5, 0)
        
    def flipY(self, pos):
        """because computer graphics has a flipped +y to the norm."""
        place = numpy.array([[1], [-1], [1]])*(pos) 
        place[1,0] += size[1]
        return place
        
    def ints(self, pos):
        """Pygame only draws on integer pixels"""
        return (int(pos[0,0]), int(pos[1,0]))
        
    def line(self, start, end):
        """draws a line from start to end, coords local to origin"""
        start = self.ints(self.flipY(start))
        end = self.ints(self.flipY(end))
        thickness = 1
        pygame.draw.aaline(self.s, self.colours["line"], start, end, thickness)
    
    def blob(self, pos):
        """Draws a blob at pos, coords local to origin"""
        pos = self.ints(self.flipY(pos))
        radius = 5
        pygame.draw.circle(self.s, self.colours["blob"], pos, radius, 0)

    def render(self):
        self.screen.blit(self.s,self.rect)
        

    



def linPath1(start, end, steps):
    x = numpy.linspace(start[0], end[0], steps)
    y = numpy.linspace(start[1], end[1], steps)
    z = numpy.linspace(0, 0, steps)
    targets = []
    for i in xrange(steps):
        targets.append((x[i], y[i], z[i]))
    return targets
    

class VirtualScreen():
    """A fake display that you write to as if it were a normal screen, but it is not
    Use for showing more than one pygame app at a time, without massive recoding
    """
    def __init__(self, screen, size, loc, gloSize):
        self.size = size
        self.loc = loc
        self.screen = screen
        self.globalSize = gloSize
        
        self.s = pygame.Surface(self.globalSize, flags=0)
        self.rect = self.s.get_rect()
        
    
    def eventHandle(self):
        pass
    
    def scale(self):
        scaled = pygame.transform.scale(self.s, self.size)
        rect = scaled.get_rect()
        rect.topleft = self.loc
        return (scaled, rect)
    
    def translate(self):
        pass
    
    def blit(self, scren, rect):
        self.s.blit(scren, rect)
        
    def render(self):
        s, r = self.scale()
        self.screen.blit(s, r)
        
    def get_size(self):
        return self.globalSize


#Lets BEGIN :D
def begin():
    
    siz = (size[0]/2, size[1]/2)
    vscreens = []
    arms = []
    x = 0
    y = 0
    for i in xrange(3):
        vscreens.append(VirtualScreen(screen, siz, (x, y), size))
        x += siz[0]
        if x >= size[0]:
            x = 0
            y += siz[1]
        arms.append(Arm("ccd", vscreens[-1]))
    arms[0].method = "ccd"
    arms[1].method = "ikT"
    arms[2].method = "ikPinv"

    myfont = pygame.font.Font("res/Kelvinized.ttf",18)
    
    start = None
    end = None    
    
    while 1:      
        for event in pygame.event.get():
            if event.type == pygame.QUIT: sys.exit()
        screen.blit(bg,bgrect)
        

        scoreline = "DISTANCE: "
        scoreboard = myfont.render(scoreline,1,textcolor)
        screen.blit(scoreboard,scoreboard.get_rect())


        #User Input
        for event in pygame.event.get():
            if event.type == pygame.QUIT: sys.exit()

            if event.type == pygame.KEYDOWN:
                
                if event.key == pygame.K_UP:
                    for arm in arms:
                        arm.sections[1]["theta"] = arm.sections[1]["theta"] + math.pi/64
                        arm.drawArm()
                elif event.key == pygame.K_DOWN:
                    for arm in arms:
                        arm.sections[1]["theta"] = arm.sections[1]["theta"] - math.pi/64
                        arm.drawArm()
                elif event.key == pygame.K_RIGHT:
                    for arm in arms:
                        arm.sections[0]["theta"] = arm.sections[0]["theta"] - math.pi/64
                        arm.drawArm()
                elif event.key == pygame.K_LEFT:
                    for arm in arms:
                        arm.sections[0]["theta"] = arm.sections[0]["theta"] + math.pi/64
                        arm.drawArm()
                elif event.key == pygame.K_SPACE:
                    pass
                elif event.key == pygame.K_e:
                    for arm in arms:
                        arm.sections[1]["length"] += 1
                        arm.drawArm()
                elif event.key == pygame.K_q:
                    for arm in arms:
                        arm.sections[1]["length"] -= 1
                        arm.drawArm()
                    
            if event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:
                    #left click
                    print event.pos
                    for arm in arms:
                        arm.target = event.pos
                        arm.inverseK(event.pos)
                    """ moving line
                    if not start:
                        start = event.pos
                    else:
                        end = event.pos
                        linpath = linPath1(start, end, 50)
                        start = None
                        for pos in linpath:
                            arm.target = (int(pos[0]), int(pos[1]))
                            arm.ik(pos[:2])
                            arm.render()
                            pygame.display.flip()
                    """
                elif event.button == 3:
                    #right click
                    for arm in arms:
                        arm.target = event.pos
                        arm.ik(event.pos)
                elif event.button == 2:
                    #middle click
                    for arm in arms:
                        arm.target = event.pos
                        arm.ccd(event.pos)#, thetaStep=5*math.pi/180)
               
                    
        for arm in arms:
            arm.render()
        for scren in vscreens:
            scren.render()

    
        pygame.display.flip()

def main():

    #draw the welcome screen
    welcome = pygame.image.load("res/welcome.png")
    wrect = welcome.get_rect()
    wrect = wrect.move(40,40)

    #wait till the user presses "enter" key
    while 1:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: sys.exit()
    
        screen.fill((55,200,44))
        screen.blit(welcome,wrect)
        pressed = pygame.key.get_pressed()
        
        if pressed[pygame.K_RETURN]: break
    
        pygame.display.flip()
    
    #BEGIN THE GAME :D
    begin()

if __name__ == "__main__":
    main()
