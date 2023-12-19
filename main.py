
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import time
import math
import numpy as np
import numpy.linalg as npla
import scipy
from scipy import linalg as spla
import scipy.sparse
import scipy.sparse.linalg
from scipy import integrate
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import warnings
np.set_printoptions(precision = 4)

#I confirm that I did not use codes from anyone else and that the work I submit is my own and my own only

def updatePos(p, v, dt):
    xNew = p[0] + (dt * v[0])
    yNew = p[1] + (dt * v[1])

    pNew = np.array([xNew,yNew])
    return pNew

def checkWallColl(p, RW, LW, TW, LoW, radius):
    if p[0]-radius > LW and p[0]+radius < RW and p[1]-radius > LoW and p[1]+radius < TW:
        return False
    else:
        return True
    
def wallCollDt(x, y, xOld, yOld, RW, LW, LoW, TW, u, v, radius):   
    dtNew = 0
    type = ""
    #epsilon = 0

    if x-radius < LW:
        #dtNew = (LW - (xOld + radius)) / u
        dtNew = (LW - (xOld - radius)) / u
        type = "LW"

    if x+radius > RW:
        dtNew = (RW - (xOld + radius)) / u
        type = "RW"

    if y-radius < LoW:
        #dtNew = (LoW - (yOld + radius)) / v
        dtNew = (LoW - (yOld - radius)) / v
        type = "LoW"

    if y+radius > TW:
        dtNew = (TW - (yOld + radius)) / v
        type = "TW"
    
    return dtNew, type

def updateVelWall(alpha, beta, v, type):
    uNew = v[0]
    vNew = v[1]
    if type == "RW" or type == "LW":
        uNew = -alpha * v[0]
        vNew = beta * v[1]
    if type == "TW" or type == "LoW":
        uNew = beta * v[0]
        vNew = -alpha * v[1]
    
    
    velNew = np.array([uNew,vNew])

    return velNew

def checkBallColl(xR, yR, xB, yB, radius):
    if math.sqrt((xB-xR)**2 +(yB-yR)**2) < (2*radius):
        return True
    else:
        return False
    
def ballCollDt(xR, yR, xB, yB, radius, uR, vR, uB, vB):
    dtShrink = abs(((math.sqrt((xB-xR)**2 + (yB-yR)**2)) - (2 * radius)) / (math.sqrt((uB-uR)**2 + (vB-vR)**2)))
    print(uR)
    print(vR)
    print(uB)
    print(vB)
    print()
    print(dtShrink)
    print()
    return dtShrink

def updateVelBall(xR, yR, xB, yB, velR, velB):
    n = np.array([xB-xR,yB-yR])
    n /= np.linalg.norm(n)
    #n1 = (xB - xR) / (math.sqrt((xB-xR)**2 + (yB-yR)**2))
    #n2 = (yB - yR) / (math.sqrt((xB-xR)**2 + (yB-yR)**2))
    #n = np.array([n1,n2])
    t = np.array([n[1], -n[0]])
    
    #uR = np.dot(velB,n)
    #vR = np.dot(velR,t)
    #uB = np.dot(velR,n)
    #vB = np.dot(velB,t)
    #vROut = np.array([uR,vR])
    #vBOut = np.array([uB,vB])
    vROut = np.dot(velB,n)*n + np.dot(velR,t)*t
    vBOut = np.dot(velR,n)*n + np.dot(velB,t)*t
    
    return vROut, vBOut

def animate(i):
    t = 0
    t_final = 10
    dt = .02
    radius = .05
    alpha = .8
    beta = .98
    vR = np.array([-.1,.5])
    vB = np.array([.11,.2])
    pR = np.array([.75,5*radius])
    pB = np.array([.25,5.5*radius])
    pROld = np.array([.75,5*radius])
    pBOld = np.array([.25,5.5*radius])
    LW = 0
    RW = 1
    TW = 1
    LoW = 0
    #xR_animation.append(pR[0])
    #yR_animation.append(pR[1])
    #xB_animation.append(pB[0])
    #yB_animation.append(pB[1])
    while t < t_final:
        t = t + dt
        #xB_animation.append(xBOld)
        #yB_animation.append(yBOld)

        checkRW = checkWallColl(pR,RW,LW,TW,LoW,radius)
        checkBW = checkWallColl(pB,RW,LW,TW,LoW,radius)
        check = checkBallColl(pR[0],pR[1],pB[0],pB[1],radius)
        
        if checkRW == True:
            dtNew, type = wallCollDt(pR[0],pR[1],pROld[0],pROld[1],RW,LW,LoW,TW,vR[0],vR[1],radius)
            pR = updatePos(pROld,vR,dtNew)
            pB = updatePos(pBOld,vB,dtNew)
           
            vR = updateVelWall(alpha,beta,vR,type)

        if checkBW == True:
            dtNew, type = wallCollDt(pB[0],pB[1],pBOld[0],pBOld[1],RW,LW,LoW,TW,vB[0],vB[1],radius)
            pB = updatePos(pBOld,vB,dtNew)
            pR = updatePos(pROld,vR,dtNew)

            vB = updateVelWall(alpha,beta,vB,type)

        if check == True:
            dtShrink = ballCollDt(pROld[0], pROld[1], pBOld[0], pBOld[1], radius, vR[0], vR[1], vB[0], vB[1])
            pR = updatePos(pROld, vR, dtShrink)
            pB = updatePos(pBOld, vB, dtShrink)
    
            vR, vB= updateVelBall(pR[0],pR[1],pB[0],pB[1],vR,vB)
        
            
        #pR = updatePos(pR,vR,dt)
        xR_animation.append(pR[0])
        yR_animation.append(pR[1]) 
        xB_animation.append(pB[0])
        yB_animation.append(pB[1])
        
        pROld = pR
        pBOld = pB
        
        pB = updatePos(pB,vB,dt) 
        pR = updatePos(pR,vR,dt)  
        

    
    ax.clear()
    ax.set_aspect(1)
    circleR = plt.Circle((xR_animation[i], yR_animation[i]), 0.05, color="red")
    circleB = plt.Circle((xB_animation[i], yB_animation[i]), 0.05, color="blue")

    ax.add_artist(circleR)
    ax.add_artist(circleB)
    ax.set_facecolor("forestgreen")
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])



# create empty lists for the x and y coordinates
xR_animation = []
yR_animation = []
xB_animation = []
yB_animation = []


# create the figure and axes objects
fig, ax = plt.subplots()
# run the animation
ani = FuncAnimation(fig, animate, frames=300, interval=100, repeat=False)

plt.show()

RW = 1
LW = 0
LoW = 0
TW = 1
for i in xR_animation:
    if i > RW or i < LW:
        print("FAIL RW or LW")
for i in yR_animation:
    if i > TW or i < LoW:
        print("FAIL TW or LoW")

for i in range(len(xR_animation)):
    if math.sqrt((xB_animation[i]-xR_animation[i])**2 +(yB_animation[i]-yR_animation[i])**2) < (2*.05):
        print("FAIL BALL COLLISION")

    #i = i + 1

