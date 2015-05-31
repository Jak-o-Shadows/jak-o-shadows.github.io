# -*- coding: utf-8 -*-
"""
Created on Thu Jul 03 17:38:11 2014

@author: Jak
"""

import math
import csv

import pygame
pygame.init()


import numpy

                 
title = "Tomodachi Life"



"""BELOW HERE IS CODE"""

def readCSVRelationships(fname):
    with open(fname, "rb") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        names = reader.next()[1:]
        relationships = []
        for i in xrange(len(names)):
            row = reader.next()[1:]
            relationships.append([int(x) for x in row])
        return names, relationships



colours = {'blue': (30, 144, 255, 255), 'brown': (58, 41, 0, 255), 'nearBlack': (26, 26, 26, 255), 'gold': (139, 105, 20, 255), 'darkblue': (0, 0, 255, 255), 'aqua': (40, 213, 146, 255), 'seagreen': (0, 72, 26, 255), 'greeny': (144, 238, 144, 255), 'grey': (127, 127, 127, 255), 'yellow': (255, 255, 0, 255), 'lightGreen': (4, 255, 0, 255), 'purple': (160, 32, 240, 255), 'pink': (255, 192, 203, 255), 'black': (0, 0, 0, 255), 'orange': (255, 165, 0, 255), 'green': (0, 41, 6, 255), 'white': (255, 255, 255, 255), 'red': (255, 0, 0, 255)}

def makeText(colour, size, text, bgcolour, textSize=15):
    """Returns a surface of colour bgcolour with text on it.
    """
    sx = int((len(text)+1)*textSize/2.5)
    size = (sx, size[1])
    image = pygame.Surface(size)
    image.fill(bgcolour)
    font = pygame.font.SysFont(None, textSize)
    txtSurface = font.render(text, False, colour, bgcolour)
    tx = (image.get_width() - txtSurface.get_width())/2
    image.blit(txtSurface, (tx, size[1]/2))
    image.convert()
    return image

def drawLine(screen, colourStart, colourEnd, start, end, thick1, thick2):
    if colourStart == colourEnd:
        pygame.draw.line(screen, colourStart, start, end, thick1)
    else:
        steps = 10
        lineX = numpy.linspace(start[0], end[0], steps)
        lineY = numpy.linspace(start[1], end[1], steps)
        r = numpy.linspace(colourStart[0], colourEnd[0], steps)
        g = numpy.linspace(colourStart[1], colourEnd[1], steps)
        b = numpy.linspace(colourStart[2], colourEnd[2], steps)
        a = numpy.linspace(colourStart[3], colourEnd[3], steps)
        for seg in xrange(steps-1):
            s = (lineX[seg], lineY[seg])
            e = (lineX[seg+1], lineY[seg+1])
            colour = (r[seg], g[seg], b[seg], a[seg])
            pygame.draw.line(screen, colour, s, e, max(thick1, thick2))





if __name__ == "__main__":
    
    names, relationShips = readCSVRelationships("tomodachi life.csv")
    
    width, height = 1000.0, 1000.0
    panel = 300
    
    #pygame screen setup        
    screen = pygame.display.set_mode([int(width + panel), int(height)])
    mapImage = pygame.Surface((width, height))
    mapImage.fill(colours["black"])
    screen.blit(mapImage, (0, 0))
    pygame.display.flip()
    pygame.display.set_caption(title)
    
    
    # customise this big to change looks
    radius = 350
    labelOffset = 75
    lineThickness = [0, 2, 2, 2, 2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
    lineColours = [(0, 0, 0, 255),#not friends / self  #0
                   (47, 21, 112, 255),#not much fun     #1
                   (153, 57, 19, 255),#mostly fun       #2
                   (204, 142, 93, 255),#fun             #3
                   (152, 178, 116, 255),#a lot of fun   #4
                   (87, 122, 38, 255),#great fun        #5
                   (63, 87, 35, 255),#most enjoyable    #6
                   (0, 0, 150, 255),#BF trust completely    #7
                   (0, 0, 255, 255),#BFF                #8
                   (107, 0, 130, 255),#very much in love      #9
                   (159, 0, 193, 255),#Completely in love     #10    
                   (210, 0, 255, 255),#want to get married      #11
                   (150, 150, 150, 255),#married - happy
                   (245, 245, 245, 255),#married - completley happy
                   (255, 255, 255, 255)] #married - ???
                   
                   
    legendLabels = ["Not much fun",
                    "Mostly fun",
                    "Fun",
                    "A lot of fun",
                    "Great fun",
                    "Most enjoyable",
                    "Best Friends-trust completely",
                    "Best Friends Forever",
                    "Very much in love",
                    "Completely in love",   
                    "Want to get married",
                    "Married - Happy",
                    "Married - Completely Happy",
                    "Married - Next stage?"]
    legendPosx = int(width + 25)
    legendYStart = 50
    legendGap = 50
    legendXgap = 50
    legendColourLength = 50
    legendColourThickness = 10
    
    
    numPeople = len(names)
    angularDistance = math.radians(360/numPeople)
    angles = [x*angularDistance for x in xrange(numPeople)]
    angles = angles[::-1]
    positions = [(int(radius*math.cos(x) +width/2), int(radius*math.sin(x)+height/2)) for x in angles]
    labelRadius = radius + labelOffset
    labelPositions = [(int(labelRadius*math.cos(x) +width/2), int(labelRadius*math.sin(x)+height/2)) for x in angles]
    
    
    
    markerRadius = 10
    
    for i in xrange(numPeople):
        text = makeText(colours["black"], (50, 20), names[i], colours["white"])
        screen.blit(text, labelPositions[i])
        for j in xrange(i, numPeople):
            if j != i:
                start = positions[i]
                end = positions[j]
                key = relationShips[i][j]
                key2 = relationShips[j][i]
                drawLine(screen, lineColours[key], lineColours[key2], start, end, lineThickness[key], lineThickness[key2])
        pygame.draw.circle(screen, colours["white"], positions[i], markerRadius)
                
                
    for i in xrange(len(legendLabels)):
        y = legendYStart + i*legendGap
        text = makeText(colours["black"], (150, 20), legendLabels[i], colours["white"])
        screen.blit(text, (legendPosx, y))
        lineStart = (legendPosx + 150 + legendXgap, y)
        lineEnd = (legendPosx + 150 + legendXgap + legendColourLength, y)
        pygame.draw.line(screen, lineColours[i+1], lineStart, lineEnd, legendColourThickness)
    
    titleText = makeText(colours["white"], (150, 100), title, (0, 0, 0, 0), textSize=75)
    screen.blit(titleText, (-25, -25))
    
    
    pygame.image.save(screen, title+".png")
    keepGoing = False
    clock = pygame.time.Clock()
    while keepGoing:
        clock.tick(30)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                keepGoing = False
                
                
        pygame.display.update()