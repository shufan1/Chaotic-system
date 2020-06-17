# Sensitivity to initial conidiions for DDP
#Author: Shufan Xia
#Date: May,2020

# import necessary packag
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import odeint # import scipy ode solver

#Use Latex fonts
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

## set parameters
omegaD = 2*np.pi
omegaN = 3*np.pi
beta = omegaN/4 #damp coefficient: b/m

############################## Linear case ##########################################
gamma = 0.1 #dirve strength coefficent tau0/ml^2 = F0/ml


# initial condition #1: A
theta0_A = 0
omega0_A = 0
r0_A = [theta0_A,omega0_A]

t = np.arange(0,10*Td,0.01)
## define the coupled DEQ with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve DEQ
r_A = odeint(f,r0_A,t)

theta_A = r_A[:,0]
omega_A = r_A[:,1]

# initial condition #2: B
theta0_B = 0.1 # initial condition slightly differed by 0.1 rad
omega0_B = 0
r0_B = [theta0_B,omega0_B]

t = np.arange(0,10*Td,0.01)
## define coupled DEQ with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve DEQ
r_B = odeint(f,r0_B,t)
theta_B = r_B[:,0]
omega_B = r_B[:,1]

delta_theta = np.abs(theta_B-theta_A) # difference between thetaB-thetaA

# plot result
plt.plot(t[:750],np.log(delta_theta[:750]))# plot in log scale
### add an approximate linear line that fit the successive maxima
x = np.arange(0,8,0.5)
y = -26/11*x-1.8
plt.plot(x,y,'--')
plt.xlabel('$t$')
plt.ylabel('log$|\Delta \\theta(t)|$')
plt.savefig("0.1deltaPhi.pdf")
plt.show()

############################# Nonlinear case ############################################################
gamma = 1.105 #dirve strength coefficent tau0/ml^2 = F0/ml


# initial condition #1: A
theta0_A = 0
omega0_A = 0
r0_A = [theta0_A,omega0_A]

t = np.arange(0,20*Td,0.01)

## define coupled DEQ with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve ODE 
r_A = odeint(f,r0_A,t)
theta_A = r_A[:,0]
omega_A = r_A[:,1]

# initial condition #2: B
theta0_B = 0.00001 # initial conditions differ even smaller, just a tiny bit
omega0_B = 0
r0_B = [theta0_B,omega0_B]

t = np.arange(0,20*Td,0.01)
## define coupled DEQ with new parameters
def f(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma*omegaN**2*np.cos(omegaD*t)

    return [ftheta,fomega]
# solve ODE 
r_B = odeint(f,r0_B,t)
theta_B = r_B[:,0]
omega_B = r_B[:,1]


delta_theta = np.abs(theta_B-theta_A) # difference between thetaB-thetaA

# plot the result
plt.plot(t[:1600],delta_theta[:1600])
plt.yscale("log")# plot in log scale
plt.axhline(1,ls='--',linewidth=0.7)
plt.xlabel('$t$')
plt.ylabel('log$|\Delta \\theta(t)|$')
plt.savefig("1.105deltaPhi.pdf")



