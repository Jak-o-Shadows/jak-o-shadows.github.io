#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pygame
import time

import matplotlib.pyplot as plt

class Joystick:
    def __init__(self, jNo):
        #pygame.joystick.init()
        self.joystick = pygame.joystick.Joystick(jNo)
        self.joystick.init()
        
        self.time = time.time()
        self.setup()
        
    def setup(self):
        self.numAxis = self.joystick.get_numaxes()
        self.data = [[] for x in xrange(self.numAxis)]
        

    def log(self):
        tim = time.time()
        for i in xrange(self.numAxis):
            axis = self.joystick.get_axis(i)
            self.data[i].append([axis, tim])

        
    def returnData(self):
        return (self.data, self.time)

def plot(data):
    d, tim = data
    dat = []
    for i in d:
        v = [j[0] for j in i]
        t = [j[1] for j in i]
        t = [j - tim for j in t]
        dat.append([v, t])
        
    for i in xrange(len(dat)):
        plt.figure(i+1)
        j = dat[i]
        plt.plot(j[1], j[0])
        plt.ylabel("Axis Value")
        plt.xlabel("Time (s)")
        plt.title("Axis " + str(i))
    plt.show()

def main():
    pygame.joystick.init()
    jTot = str(pygame.joystick.get_count())
    print "Starting from zero enter the joystick number. There are "+jTot+" Joysticks"
    jNo = int(raw_input(":"))
    clock = pygame.time.Clock()
    FPS = 30
    joystick = Joystick(jNo)
    
    screen = pygame.display.set_mode((100, 100))
    
    keepGoing = True
    while keepGoing:
        clock.tick(FPS)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                keepGoing = False
        joystick.log()
    plot(joystick.returnData())
    

if __name__ == '__main__':
    main()