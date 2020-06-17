
# Bifurcation of DDP
#Author: Shufan Xia
#Date: May,2020

# import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import odeint # import scipy ode solver

#Use Latex fonts
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

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

# ---------------------------------------------------------------------------------------------------

# set parameters
omegaD =2*np.pi
omegaN = 3*np.pi
beta = omegaN/4 #damp coefficient: b/m
Td= 2*np.pi/omegaD

# x-axis of the bifurcation diagram. gamma values range from 1.06 to 1.087
gamma = np.linspace(1.06,1.087,271,endpoint=True) # 

#same initial conditions p500 Taylor
theta0=-np.pi/2
omega0=0 
r0 = [theta0,omega0]

critical_gamma = []
periodicity_0 = 1
first_time = False
# collect data for x and y values to plot bifurcation diagram
x=[]
y=[]

for gamma_i in gamma:

    t = np.arange(0,700*Td,0.01) ## from t=0 to 700s

    # define DEQ
    def f(r,t):
        theta = r[0]
        omega = r[1]
        ftheta = omega
        fomega = -omegaN**2*np.sin(theta)-2*beta*omega+gamma_i*omegaN**2*np.cos(omegaD*t)
        return [ftheta,fomega]
    # solve DEQ
    r = odeint(f,r0,t)
    theta = r[:,0]
    omega = r[:,1]
    
    theta=correctTheta(theta[50100:55010][::100])
    
    # set preceision to 0.001
    theta = np.round(theta,3)
    
    
    theta_0 =  theta[0]
    theta_t = [theta_0]
    periodicity = 1
    
    for i in range(len(theta)):
        
        # count number of differnt theya_t to find periodicity
        if not theta[i] in theta_t:
            periodicity +=1
            theta_t.append(theta[i])        
        ### add to the data will be used to plot bifurcation diagram
        x.append(gamma_i)
        y.append(theta[i])
        
    # pick up threshold value of bifurcation
    if ((periodicity==2) or (periodicity==4) or (periodicity==8) or (periodicity==15) )and (periodicity != periodicity_0):
        periodicity_0 = periodicity
        critical_gamma.append(gamma_i)
        
        
        
# Plot bifucration diagram
f,ax = plt.subplots(1)
ax.plot(x,y,'.',markersize=1)
ax.axvline(critical_gamma[0],ymin=0,ymax=0.74,ls='--',linewidth='1.5')
ax.axvline(critical_gamma[1],ymin=0,ymax=0.9,ls='--',lw='1.5')
ax.axvline(critical_gamma[2],ymin=0,ymax=0.95,ls='--',lw='1.5')
ax.axvline(critical_gamma[3],ymin=0,ymax=0.95,ls='--',lw='1.5')
# add another x-axis label for critical gamma 
ax2 = ax.twiny()
ax2.set_xticks([1.0662,1.0794,1.0821,1.0831])
ax2.set_xticklabels(['$\gamma_1$', '$\gamma_2$', '$\gamma_3$','$\gamma_4$'])
ax2.xaxis.set_ticks_position('bottom')
ax2.set_xlim(ax.get_xlim())
# adjust the original x-axis
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_coords(0.5,-0.1)
ax.set_xlabel('$\gamma$')
ax.set_ylabel('$\\theta(t)$')

plt.savefig("bifurcate.pdf")
plt.show()



# Find Feigenbaum number
critical_gamma=[1.0663,1.0793,1.0822,1.0828]
interval = [""]
for i in range(len(critical_gamma)-1):
    interval.append(np.round(critical_gamma[i+1]-critical_gamma[i],4))

delta1=interval[1]/interval[2]
delta2=interval[2]/interval[3]
print(interval)
print(delta1,delta2)
print("Feigenbaum  number is: ",(delta1+delta2)/2)