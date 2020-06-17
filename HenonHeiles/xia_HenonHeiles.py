# Henon-Heiles system: solution, 3D state orbit
#Author: Shufan Xia
#Date: May,2020

import numpy as np
# import numpy and matplot library
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
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


###############################################################################################################
#### plot the countour plot of HH potential
x,y=np.meshgrid(np.arange(-1,1.1,0.1),np.arange(-1,1.4,0.1))
potential = V(x,y)

fig,ax=plt.subplots(1,1)
cp = ax.contour(x, y, potential,[1/100,1/24,1/12,1/8,1/6])
ax.clabel(cp, fontsize=12)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.savefig("potential.pdf")
plt.show()


###############################################################################################################
############################################## E = 1/12 #######################################################
# E=1/12, x0=0,y0=0, ydot0=0,y xdot0??
x0 = 0
y0 = 0
ydot0 = 0
E = 1/12
xdot0 = fsolve(solve_xdot0, 0.5) # solve xdot0 that satisfy the inital condition

r0 = [x0,xdot0,y0,ydot0]
t = np.arange(0,1000,0.01)

# solve DEQ
r = odeint(f,r0,t)
x = r[:,0]
xdot = r[:,1]
y = r[:,2]
ydot = r[:,3]

#### plot solution to x and y
plt.subplot(211)
plt.plot(t[:16000],x[:16000])
plt.xlabel("$t$")
plt.ylabel("$x$")
plt.subplot(212)
plt.plot(t[:16000],y[:16000])
plt.xlabel("$t$")
plt.ylabel("$y$")
plt.savefig("1-12Sol.pdf")
plt.show()

# plot 3D phase plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x[:50000], y[:50000], ydot[:50000], linewidth=0.4)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$\dot{y}$")
plt.savefig("1-123D.pdf")
plt.show()

# plot Poincare section
# two arrary to hold x and y value for plotting
yP = []
ydotP=[]

# select y and ydot at which x=0 and x
for i in range(len(x)-1):
    if (x[i]<=0 and x[i+1] >=0):
    	# take the average of y and ydot in the neighboring values 
        yP.append((y[i]+y[i+1])/2) 
        ydotP.append((ydot[i]+ydot[i+1])/2)

plt.plot(yP,ydotP,'.',markersize=3,label=ydot0)
plt.xlabel("$y$")
plt.ylabel("$\dot{y}$")
plt.savefig("1-12Poin_one.pdf")
plt.show()

###############################################################################################################
############################################## E = 1/8 #######################################################
#E=1/8, x0=0,y0=0, ydot0=0,y xdot0??
x0 = 0
y0 = 0
ydot0 = 0
E = 1/8
xdot0 = fsolve(solve_xdot0, 0.5) # solve xdot0 that satisfy the inital condition

r0 = [x0,xdot0,y0,ydot0]
t = np.arange(0,1000,0.01)

# solve DEQ
r = odeint(f,r0,t)
x = r[:,0]
xdot = r[:,1]
y = r[:,2]
ydot = r[:,3]

# plot 3D phase plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x[:50000], y[:50000], ydot[:50000], linewidth=0.4)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$\dot{y}$")
plt.savefig("1-83D.pdf")
plt.show()

###############################################################################################################
############################################## E = 1/6 #######################################################
# E=1/6, x0=0,y0=0, ydot0=0,y xdot0??
x0 = 0
y0 = 0
ydot0 = 0
E = 1/6
xdot0 = fsolve(solve_xdot0, 0.5) # solve xdot0 that satisfy the inital condition

r0 = [x0,xdot0,y0,ydot0]
t = np.arange(0,1000,0.01)

# solve DEQ
r = odeint(f,r0,t)
x = r[:,0]
xdot = r[:,1]
y = r[:,2]
ydot = r[:,3]

# plot 3D phase plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x[:50000], y[:50000], ydot[:50000], linewidth=0.4)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$\dot{y}$")
plt.savefig("1-63D.pdf")
plt.show()