#Poincare section of Henon-Heiles system
#Author: Shufan Xia
#Date: May,2020

import numpy as np
# import numpy and matplot library
import matplotlib.pyplot as plt
# import scipy ode solver and fsolve
from scipy.integrate import odeint
from scipy.optimize import fsolve

#Use Latex fonts
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

lam = 1 # https://jfuchs.hotell.kau.se/kurs/amek/prst/11_hehe.pdf

# define Henon-Heiles potential
def V(x,y): 
    v = 1/2*(x**2+y**2)+lam*(x**2*y-y**3/3)
    return v

# total energy of the system as a c
def E(x0,y0,xdot0,ydot0):
    e=V(x0,y0)+1/2*(xdot0**2+ydot0**2)
    return e

# a series of function to solve fot the initial values of x,y,xdot,ydot 
#          with total energy E.
#          Given E and the inital values for three variables, find the initial value of the other variable.

def solve_x0(x0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E
    return equation

def solve_y0(y0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E  
    return equation

def solve_xdot0(xdot0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E
    return equation

def solve_ydot0(ydot0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E
    return equation

### define the system of deqs with new parameters
def f(r,t):
    x = r[0]
    xdot = r[1]
    y = r[2]
    ydot = r[3]
    fx = xdot # first derivative of x
    fy = ydot # first derivative of y
    fxdot = - 1*(x + 2*lam*x*y) # ddx
    fydot = -1*(y+lam*(x**2 - y**2)) # ddy
   
    return [fx,fxdot,fy,fydot]


#########################################   E=1/12  ####################################################
f,ax=plt.subplots(1,figsize=(6,4))
E = 1/12

# case1: E=1/12, x=0 0,y0=0,ydot>0, vx=?
x0 = 0
y0 = 0
vy0 = np.arange(0,0.45,0.1)


for ydot0 in vy0:   
    vx0 = fsolve(solve_xdot0, [0.5,-0.5])  # find vx_0 for each vy_0
    
    yP = [] # hold data for plotting
    ydotP=[]
    
    for xdot0 in vx0:
        
        r0 = [x0,xdot0,y0,ydot0]

        t = np.arange(0,1000,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 

        for i in range(len(x)-1):
            if (x[i]<=0 and x[i+1] >=0):
                # take the average of y and ydot in the neighboring values 
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=ydot0) # plot data

# case1: E=1/12, xdot=0, ydot=0, y0>0, x0=??
xdot0 = 0
ydot0 = 0
y_initial = [0.09,0.39,0.4]

for y0 in y_initial:
    x0 = fsolve(solve_x0,[0.5,-0.5]) # solve x_0 for each y_0

    yP = [] # hold data for plotting
    ydotP = []


    for i in range(len(x0)):
        r0 = [x0[i],xdot0,y0,ydot0]

        t = np.arange(0,2500,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 
        for i in range(len(x)-1):
            if (x[i]<=0 and x[i+1] >=0):
                # take the average of y and ydot in the neighboring values 
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=y0)   # plot data

ax.set_xlabel('$y$')
ax.set_ylabel('$\dot{y}$')
plt.savefig("1-12Poin.pdf")
plt.show()

#######################################################################################################
#########################################   E=1/8  ####################################################
f,ax=plt.subplots(1,figsize=(6,4))
E = 1/8

# case1: E=1/8, x=0 0,y0=0,ydot=>0, vx=?
x0 = 0
y0 = 0
vy0 = [0.0,0.1,0.2,0.4,0.5]


for ydot0 in vy0:   
    vx0 = fsolve(solve_xdot0, [0.5,-0.5])  # find vx_0 for each vy_0
    
    yP = [] # hold data for plotting
    ydotP=[]
    
    for xdot0 in vx0:
        
        r0 = [x0,xdot0,y0,ydot0]

        t = np.arange(0,1000,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 

        for i in range(len(x)-1):
            if (x[i]<=0 and x[i+1] >=0):
                # take the average of y and ydot in the neighboring values 
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=ydot0) # plot data

# case2: E=1/8, x_0=0, ydot=0, y0>0 and xdot_0=??
x_0 = 0
ydot0 = 0
y_initial = [0.09,0.39]

for y0 in y_initial:
    xdot0 = fsolve(solve_xdot0,[0.5,-0.5]) # solve xdot_0 for each y_0

    yP = [] # hold data for plotting
    ydotP = []


    for i in range(len(xdot0)):
        r0 = [x0,xdot0[i],y0,ydot0]

        t = np.arange(0,2500,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 
        for i in range(len(x)-1):
            if (x[i]<=0 and x[i+1] >=0):
                # take the average of y and ydot in the neighboring values 
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=y0)   # plot data

ax.set_xlabel('$y$')
ax.set_ylabel('$\dot{y}$')
plt.savefig("1-8Poin.pdf")
plt.show()


#######################################################################################################
#########################################   E=1/6  ####################################################
f,ax=plt.subplots(1,figsize=(6,4))
E = 1/6

# case1: E=1/6, x=0 0,y0=0,ydot=>0, vx=?
x0 = 0
y0 = 0
vy0 = np.arange(0,0.51,0.25)


for ydot0 in vy0:   
    vx0 = fsolve(solve_xdot0, [0.5,-0.5])  # find vx_0 for each vy_0
    
    yP = [] # hold data for plotting
    ydotP=[]
    
    for xdot0 in vx0:
        
        r0 = [x0,xdot0,y0,ydot0]

        t = np.arange(0,1000,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 

        for i in range(len(x)-1):
            # take the average of y and ydot in the neighboring values 
            if (x[i]<=0 and x[i+1] >=0):
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=ydot0) # plot data

# case2: E=1/6, x_0=0, ydot_0=0, y0>0 and xdot_0=??
x_0 = 0
ydot0 = 0
y_initial = [0.2,0.4,0.5]

for y0 in y_initial:
    xdot0 = fsolve(solve_xdot0,[0.5,-0.5]) # solve xdot_0 for each y_0

    yP = [] # hold data for plotting
    ydotP = []


    for i in range(len(xdot0)):
        r0 = [x0,xdot0[i],y0,ydot0]

        t = np.arange(0,2500,0.01)

        # define the system of deqs with new parameters
        def f(r,t):
            x = r[0]
            xdot = r[1]
            y = r[2]
            ydot = r[3]
            fx = xdot                           #dx
            fy = ydot                           #dy
            fxdot = - 1*(x + 2*lam*x*y)         #ddx
            fydot = -1*(y+lam*(x**2 - y**2))    #ddu
            return [fx,fxdot,fy,fydot]

        r = odeint(f,r0,t) # solve DEQ

        x = r[:,0]
        xdot = r[:,1]
        y = r[:,2]
        ydot = r[:,3]

        # select y and ydot with x=0 and xdot > 0 
        for i in range(len(x)-1):
            if (x[i]<=0 and x[i+1] >=0):
                # take the average of y and ydot in the neighboring values 
                yP.append((y[i]+y[i+1])/2)
                ydotP.append((ydot[i]+ydot[i+1])/2)

    ax.plot(yP,ydotP,'.',markersize=1.5,label=y0)   # plot data

ax.set_xlabel('$y$')
ax.set_ylabel('$\dot{y}$')
plt.savefig("1-6Poin.pdf")
plt.show()