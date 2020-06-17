#Damped Driven Pendulum system: solution, phase plot and poincare section
#Author: Shufan Xia
#Date: May,2020


import numpy as np
# import numpy and matplot library
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits import mplot3d
# import scipy ode solver
from scipy.integrate import odeint

#Use Latex fonts
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

### case 1: small driving strength
## approximately linear


gamma = 0.2 #dirve strength coefficent tau0/ml^2 = F0/ml
omegaD = 2*np.pi  # condition in Taylor p480 w0=2*pi, so that Td=1

omegaN = 1.5*omegaD
beta =  omegaN/4 #damp coefficient: b/m, beta/2 = omegaN/4 in Taylor

# for sqrt(g/l) = 1.5 wd
g = 9.8
l = g/omegaN**2

# find driving and natural period
T0 = 2*np.pi/omegaN
Td =2*np.pi/omegaD

# set initial condition
theta0 = 0 #0.3*np.pi/2 #start with small angle
omega0 = 0
r0 = [theta0,omega0]

# 
t = np.arange(0,8*Td,0.01)

# solve DEQ
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]

r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]
drivePhase = omegaD*t

# plot solution
f=plt.figure(figsize=(13,2))
plt.plot(t,theta)
plt.xlabel("$t$")
plt.ylabel("$\\theta$")
plt.grid()
plt.savefig("0.2soultion.pdf")
plt.show()

# phase plot
plt.plot(theta,omega)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
plt.grid()
plt.savefig("0.2phase.pdf")
plt.show()

#### 3D phase plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(theta, omega, drivePhase, 'gray')
ax.set_xlabel('$\\theta$')
ax.set_ylabel('$\dot{\\theta}$')
ax.set_zlabel('$\\omega_d t$')
plt.savefig("0.23D.pdf")
plt.show()

# Poincare section, 100 elements per cycle
thetaTransient=theta[::100][:3] ##  theta abd omega when the system is the transient state
omegaTransient=omega[::100][:3]
thetaP=theta[::100][3:] ##  theta abd omega when the system after the transient state
omegaP=omega[::100][3:]
plt.plot(thetaTransient,omegaTransient,'b.')
plt.plot(thetaP,omegaP,'k.')
# specify min and max of theta and theta dot as the limit of the plot
plt.xlim(-0.4,0.4)
plt.ylim(-2.5,2.5)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
plt.grid()
plt.savefig("0.2Poin.pdf")
plt.show()

##### comment out if you need to produce a txt file of the solution
# data = np.array([t,theta,omega])
# data = data.T
# #here you transpose your data, so to have it in two columns
# datafile_path = "C:/Users/lenovo/Desktop/PHYS304/final/0_2.txt"
# with open(datafile_path, 'w+') as file:
# #here you open the ascii file
#     np.savetxt(file, data, fmt=["%.2f","%.5f","%.5f"])

#############################################################
### case 2: gamma = 0.9
### same damping coefficient, little bit larger amplitude of driving force nonlinearity shows up

gamma = 0.9 #dirve strength coefficent tau0/ml^2 = F0/ml
omegaD = 2*np.pi  # condition in Taylor p480 w0=2*pi, so that Td=1
omegaN = 1.5*omegaD
beta = omegaN/4 #damp coefficient: b/m

Td =2*np.pi/omegaD

# set initial condition
theta0 = 0 
omega0 = 0
r0 = [theta0,omega0]


t = np.arange(0,10*Td,0.01)
# define DEQ with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]

# solve ODE 
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]
drivePhase =  omegaD*t # drive phase

# plot solution
f = plt.figure(figsize=(13,2))
plt.plot(t,theta,label="$\\theta(t)$")
plt.legend()
plt.grid()
plt.xlabel("$t$")
plt.ylabel("$\\theta(t)$")
plt.savefig("0.9solution.pdf")
plt.show()

# plot phase plot
plt.plot(theta,omega)
plt.xlabel("$\\theta(t)$")
plt.ylabel("$\dot{\\theta}(t)$")
plt.savefig("0.9phase.pdf")
plt.show()

# 3D phase plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(theta, omega, drivePhase, 'gray')
ax.set_xlabel('$\\theta$')
ax.set_ylabel('$\dot{\\theta}$')
ax.set_zlabel('$\\omega_d t$')
plt.savefig("0.93D.pdf")
plt.show()

#### the soultion is not a pure sin from t=4.83 to 6.25, zoom in and fit a pure sine function
thetaPi= list(map(lambda thetai: thetai, theta[493:620])) 
# fit a pure sine function
amplitude = max(theta[483:625])
phase = np.arcsin(theta[500]/amplitude)

# theta(t) if the solution is pure sin, 
sin = list(map(lambda ti: amplitude*np.sin(omegaD*ti+phase),t[493:620]))

plt.plot(t[493:620],thetaPi,label="numerical solution")
plt.plot(t[493:620],sin,'--',label="pure sin function")
plt.xlabel("$t$")
plt.ylabel("$\\theta$")
plt.legend()
plt.grid()
plt.savefig("0.9zoom.pdf")
plt.show()

# Poincare section,
thetaTransient=theta[::100][:5]
omegaTransient=omega[::100][:5]
thetaP=theta[::100][5:]
omegaP=omega[::100][5:]
plt.plot(thetaTransient,omegaTransient,'b.') # plot transient and steady state with different colors
plt.plot(thetaP,omegaP,'k.')
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
plt.grid()
plt.savefig("0.9Poin.pdf")
plt.show()

##### comment out if you need to produce a txt file of the solution
# data = np.array([t,theta,omega])
# data = data.T
# #here you transpose your data, so to have it in two columns
# datafile_path = "C:/Users/lenovo/Desktop/PHYS304/final/0_9.txt"
# with open(datafile_path, 'w+') as file:
# #here you open the ascii file
#     np.savetxt(file, data, fmt=["%.2f","%.5f","%.5f"])
##############################################################################################

### case 3: gamma >1 =1.073

beta = omegaN/4 #damp coefficient: b/m
gamma = 1.073 #dirve strength coefficent tau0/ml^2 = F0/ml


# set initial condition
theta0 = 0#.2*np.pi/2 #start with small angle
omega0 = 0
r0 = [theta0,omega0]

# over 330 period
t = np.arange(0,330*Td,0.01)

# define DEQ
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve DEQ
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]


### zoom in after the transient state gone：to see period doubling
t_Zoom=t[3000:4000]
theta_Zoom=theta[3000:4000]
f,ax=plt.subplots(1,figsize=(5,3))
ax.plot(t_Zoom,theta_Zoom/np.pi)
plt.ylabel("$\\theta$")

### draw a line from the lowest truough
## set y-axis tick in terms of pi
## https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))


### show the effect of period doubling, doubleing twice, a full period = 2 cycle
y1=min(theta_Zoom)/np.pi# the value of the lowest trough, y value of the line across the bottom of the trough
index = np.argmin(theta_Zoom) 
t1 = t_Zoom[index] #t1,t2,t3 x-value where the trough occurs
t2 = t1+2
t3 = t1+4 
x_values = [t1,t2,t3]
y_values = [y1,y1,y1]

plt.plot(x_values, y_values,'o--',label="one complete cycle")
plt.legend()
plt.savefig("1.073zoom.pdf")
plt.show()

#---------- two functions to produce phase plot for the system after transient state is gone --------------------
# plot phase plot, 
#			t_stable: where you think the system is in 
#			useRad:  set the tick of x-axis (theta) as n*pi or not
def PhasePlot(t_stable,useRad):
    # phase plot after stable
    i = t_stable*100
    fig,ax = plt.subplots(1)

    if (useRad == "Rad"):
        ax.plot(theta[i:]/np.pi,omega[i:],'.',markersize="0.8")
    ### set x and y-axis ticker in terms of multiple of pi
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    else:
        ax.plot(theta[i:],omega[i:],'.',markersize="1")
        
    plt.xlabel("$\\theta$")
    plt.ylabel("$\dot{\\theta}$")
    plt.grid()


# plot PoinCare section of the phase plot after the transient state ends, the system enters the stable state
#			t_stable: where you think the system is in 
#			useRad:  set the tick of x-axis (theta) as n*pi or not
def PoinCare(t_stable,useRad):

    thetaP=theta[::100][t_stable:]
    omegaP=omega[::100][t_stable:]

    fig,ax = plt.subplots(1)
    ax.plot(thetaP,omegaP,'.',markersize='6')
    plt.xlabel("$\\theta$")
    plt.ylabel("$\dot{\\theta}$")
    ### set x and y-axis ticker in terms of multiple of pi
    if (useRad == "Rad"):
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    plt.grid()

#---------------------------------------------------------------------------------------------
# plot phase plot and poincare section for gamma - 1.073
PhasePlot(26,"Deg")
plt.savefig("1073phase.pdf")
plt.show()
PoinCare(26,"Deg")
plt.savefig("1073Poin.pdf")
plt.show()

##### comment out if you need to produce a txt file of the solution
# data = np.array([t,theta,omega])
# data = data.T
# #here you transpose your data, so to have it in two columns
# datafile_path = "C:/Users/lenovo/Desktop/PHYS304/final/1_073.txt"
# with open(datafile_path, 'w+') as file:
# #here you open the ascii file
#     np.savetxt(file, data, fmt=["%.2f","%.5f","%.5f"])
##############################################################################################


### case 4: gamma - 1.081
beta = omegaN/4 #damp coefficient: b/m
gamma = 1.081 #dirve strength coefficent tau0/ml^2 = F0/ml


# set initial condition
theta0 = -1*np.pi/2 #start with small angle
omega0 = 0
r0 = [theta0,omega0]


t = np.arange(0,330*Td,0.01)
# define DEQ with new parmaeters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve DEQ
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]

# plot solution
f,ax=plt.subplots(1,figsize=(13,2))
ax.plot(t[:1000],theta[:1000]/np.pi,label="$\\theta(t)$")
plt.axhline(0,color='k',linewidth=0.3)
plt.grid()
plt.xlabel("$t$")
plt.ylabel("$\\theta$")
plt.legend()
### set ticker in terms of multiple of pi
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
plt.show()


### zoom in after the transient state is gone： period doubling. t=16 t0 28, 
t_Zoom=t[1620:2820]
theta_Zoom=theta[1620:2820]
omega_Zoom = omega[1620:2820]
f,ax=plt.subplots(1,figsize=(5,3))
ax.plot(t_Zoom,theta_Zoom/np.pi,label="$\\theta(t)$")#,color='k',linewidth=0.8)
plt.xlabel("$t(s)$")
plt.ylabel("$\\theta$")
plt.grid()
### set ticker in terms of multiple of pi
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))


### to see the effect of period doubling twice draw a horizontal line at the lowest trough
index = np.argmin(theta_Zoom)
t1 = t_Zoom[index]
t2 = t1+4
y1=min(theta_Zoom)/np.pi
x_values = [t1, t2]
y_values = [y1, y1]
plt.plot(x_values, y_values,'o--',label="one complete cycle")
plt.legend()
plt.savefig("1081zoom.pdf")
plt.show()


# make phase plot and poincare section for
PhasePlot(16,"Deg")
plt.savefig("1081phase.pdf")
plt.show()

PoinCare(16,"Deg")
plt.savefig("1081poin.pdf")
plt.show()

###### see period doubling twice numerically
t_i=t_Zoom[index:][::100][:10]#[index:index+10]
theta_i=np.around(theta_Zoom[index:][::100][:10],5)
omega_i=np.around(omega_Zoom[index:][::100][:10],5)
cell_text=np.column_stack((t_i,theta_i,omega_i))
columns=['$t(s)$','$\\theta(t)$','$\dot{\\theta}(t)$']
f,ax=plt.subplots(1)
table = plt.table(cellText=cell_text,
            colLabels=columns,
            colWidths=[0.2] * 3,
            loc="center")
ax.axis("off")
table.scale(2, 2)
plt.savefig("1.081Table.pdf")
plt.show()
##### comment out if you need to produce a txt file of the solution
# data = np.array([t,theta,omega])
# data = data.T
# #here you transpose your data, so to have it in two columns
# datafile_path = "C:/Users/lenovo/Desktop/PHYS304/final/1_081.txt"
# with open(datafile_path, 'w+') as file:
# #here you open the ascii file
#     np.savetxt(file, data, fmt=["%.2f","%.5f","%.5f"])
##############################################################################################

### case 5: gamma =1.0826, periodicity 8
beta = omegaN/4 #damp coefficient: b/m
gamma = 1.0826 #dirve strength coefficent tau0/ml^2 = F0/ml


# set initial condition
theta0 = -1*np.pi/2
omega0 = 0
r0 = [theta0,omega0]


t = np.arange(0,330*Td,0.01)

# define f with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]

# solve DEQ
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]


# plot solution
f,ax=plt.subplots(1,figsize=(13,2))
ax.plot(t[:1500],theta[:1500]/np.pi,label="$\\theta(t)$")
plt.axhline(0,color='k',linewidth=0.3)
plt.grid()
plt.xlabel("$t$")
plt.ylabel("$\\theta$")
plt.legend()
### set ticker in terms of multiple of pi
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
plt.show()


### zoom in after the transient state is gone： see  the effect of periodicity 8
t_Zoom=t[3800:5020]
theta_Zoom=theta[3800:5020]
omega_Zoom = omega[3800:5020]
f,ax=plt.subplots(1,figsize=(6,5))
ax.plot(t_Zoom,theta_Zoom/np.pi,label="$\\theta(t)$")#,color='k',linewidth=0.8)
plt.xlabel("$t(s)$")
plt.ylabel("$\\theta$")
plt.grid()
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$')) ### set ticker in terms of multiple of pi
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))

# draw a horizontal line at the lowest trough
index = np.argmin(theta_Zoom)
t1 = t_Zoom[index]
t2 = t1+8
y1=min(theta_Zoom)/np.pi
x_values = [t1, t2]
y_values = [y1, y1]
plt.plot(x_values, y_values,'o--',label="one complete cycle")
plt.legend()
plt.savefig("1.0826sol.pdf")
plt.show()

PhasePlot(30,"Deg")
plt.savefig("1.0826Phase.pdf")
plt.show()
PoinCare(30,"Deg")
plt.savefig("1.0826Poin.pdf")
plt.show()
##################################################################################################################

### case 6: gamma = 1.4, see rolling motion by looking at the phase plot
omegaD =2*np.pi
omegaN = 3*np.pi
beta = omegaN/4 #damp coefficient: b/m
gamma = 1.4  #dirve strength coefficent tau0/ml^2 = F0/ml

Td= 2*np.pi/omegaD

# set initial condition
theta0 = 0
omega0 = 0
r0 = [theta0,omega0]
t = np.arange(0,300*Td,0.01)
# update f with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve ODE 
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]

# plot first 6 periods of the phase plot
f,ax=plt.subplots(1,figsize=(13,2))
ax.plot(theta[:601]/np.pi,omega[:601],'.',markersize=3)
ax.plot(theta[::100][:7]/np.pi,omega[::100][:7],'b.',markersize=9)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
plt.grid()
# ### set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
plt.savefig("14roll.pdf")
plt.show()


#---------------------------------------------------------------------------------------
# a function to correct theta if |theta|>pi so that theta is within -pi to pi
@ np.vectorize
def correctTheta(angle):
    if(angle>=0):
        a= angle%(2*np.pi)
        if a <=np.pi:
            return a/np.pi
        else:# t>1
            return -1*(a-np.pi)/np.pi
    else:
        a= angle % (-2*np.pi)
        if a >=-1*np.pi:
            return a/np.pi
        else: #a<-1:
            return (2*np.pi+a)/np.pi

#-------------------------------------------------------------------------------------------------------------------
### change plot style, emphasize x=0 and y=0 axies with arrows in their positive direction
def show_axis():
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # make arrows
    ax.plot((1), (0), ls="", marker=">", ms=5, color="k",
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot((0), (1), ls="", marker="^", ms=5, color="k",
            transform=ax.get_xaxis_transform(), clip_on=False)


####################################################################################################################333
## case 7: gamma = 1.5, chaos 
## big and fast rotation, need to restrain theta in -pi and pi 


omegaD =2*np.pi
omegaN = 3*np.pi
beta = omegaN/4 #damp coefficient: b/m
gamma = 1.5  #dirve strength coefficent tau0/ml^2 = F0/ml

Td= 2*np.pi/omegaD

# set initial condition
theta0 = 0
omega0 = 0
r0 = [theta0,omega0]

t = np.arange(0,5000*Td,0.01)
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve ODE 
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]

theta_Adjusted = correctTheta(theta)

# -------------phase plot for 190s----------------
f,ax=plt.subplots(1)

ax.plot(theta_Adjusted[1000:20000],omega[1000:20000],'.',markersize="0.8")

plt.axhline(0,color='k',linewidth=0.4)
plt.axvline(0,color='k',linewidth=0.4)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
# set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
show_axis() # make x and y axis more obvious

plt.savefig("1.5phase.pdf")
plt.show()


# ---------poincare section --------------

# select theta nad omega every once 100 data points (1 cycle) after transient state is gone
thetaP=theta_Adjusted[::100][10:]
omegaP=omega[::100][10:]

f,ax=plt.subplots(1,figsize=(6,4))

ax.plot(thetaP,omegaP,'.',markersize='0.8')

plt.axhline(0,color='k',linewidth=0.4)
plt.axvline(0,color='k',linewidth=0.4)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
# ### set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
show_axis()# make x and y axis more obvious
plt.savefig("1.5Poin.pdf")
plt.show()

##### comment out if you need to produce a txt file of the solution
# data = np.array([t,theta,omega])
# data = data.T
# #here you transpose your data, so to have it in two columns
# datafile_path = "C:/Users/lenovo/Desktop/PHYS304/final/1_5.txt"
# with open(datafile_path, 'w+') as file:
# #here you open the ascii file
#     np.savetxt(file, data, fmt=["%.2f","%.5f","%.5f"])


########################################################################################################################3
## case 8: gamma = 1.5 but with smaller damping chaos 
## big and fast rotation, need to restrain theta in -pi and pi 


omegaD =2*np.pi
omegaN = 3*np.pi
beta = omegaN/8 #damp coefficient: b/m
gamma = 1.5  #dirve strength coefficent tau0/ml^2 = F0/ml

Td= 2*np.pi/omegaD

# set initial condition
theta0 = 0
omega0 = 0
r0 = [theta0,omega0]

t = np.arange(0,5000*Td,0.01)
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve ODE 
r = odeint(f,r0,t)
theta = r[:,0]
omega = r[:,1]

theta_Adjusted = correctTheta(theta)

# -------------phase plot for 190s----------------
f,ax=plt.subplots(1)

ax.plot(theta_Adjusted[1000:20000],omega[1000:20000],'.',markersize="0.8")

plt.axhline(0,color='k',linewidth=0.4)
plt.axvline(0,color='k',linewidth=0.4)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
# set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
show_axis() # make x and y axis more obvious

plt.savefig("1.5dampphase.pdf")
plt.show()


# ---------poincare section --------------

# select theta nad omega every once 100 data points (1 cycle) after transient state is gone
thetaP=theta_Adjusted[::100][10:]
omegaP=omega[::100][10:]

f,ax=plt.subplots(1,figsize=(6,4))

ax.plot(thetaP,omegaP,'.',markersize='0.8')

plt.axhline(0,color='k',linewidth=0.4)
plt.axvline(0,color='k',linewidth=0.4)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
# ### set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
show_axis()# make x and y axis more obvious
plt.savefig("1.5dampPoin.pdf")
plt.show()

#-------------------zoom in at the upper right region ----------------
f,ax=plt.subplots(1,figsize=(6,4))
ax.plot(thetaP,omegaP,'.',markersize='0.8')
# specify min and max of theta and theta dot as the limit of the plot

plt.xlim(0,1)
plt.ylim(12,30)
plt.xlabel("$\\theta$")
plt.ylabel("$\dot{\\theta}$")
# set ticker in terms of multiple of pi
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
plt.savefig("1.5dampZoom.pdf")
plt.show()
