# author: Shufan Xia
# date: May,2020

from vpython import *
# import numpy and matplot library
import matplotlib.pyplot as plt
import numpy as np

filename=input("filename:")
data=np.loadtxt(filename,float)

# read in data file
t = data[:,0]
theta = data[:,1]
theta0=theta[0]
omega = data[:,2]


l=0.1

x0 = l*cos(theta0-np.pi/2)
y0 = l*sin(theta0-np.pi/2)+0.07

## create 3 plot windows

# 1 theta(t), omega(t)
graph1 = graph(xtitle='time', ytitle='value', scroll=True, fast=False, xmin=0,xmax=8, width=1100, height=200)
funct1 = gcurve(color=color.blue, width=4, label='theta')
funct2 = gcurve(color=color.orange, width=4, label='omega')

# 3 poincare
graph3 = graph(x=550, y=0, title="Poincare", xtitle='time', ytitle='value', fast=False, width=450,height=360,align='right')
funct4 = gdots(color=color.blue, radius=2, label='dots')

# 2 phase
graph2 = graph(x=0,y=0,title="phase plot", xtitle='time', ytitle='value', fast=False, width=450,height=360,align='left')
funct3 = gdots(color=color.red, radius=0.8, label='dots')

# initial position for the sphere


for i in range(len(theta)):

	if i< 7000:

		rate(500)
		th = theta[i]
		dth = omega[i]*0.01

		x = l*cos(th-np.pi/2)
		y = l*sin(th-np.pi/2)+0.07

		# update the position of the vector
		# v.rotate(dth,axis=vector(0,0,1),origin=vector(0,0.07,0))
		# # update the position of the shphere
		# s.pos=vector(x,y,0)

		funct1.plot(t[i], theta[i])
		funct2.plot(t[i], omega[i])

		if i >= 1200:
			funct3.plot(theta[i],omega[i])

			if i%100 ==0:
				funct4.plot(theta[i],omega[i])

	if i >= 7000 and i <=50000:
		rate (1500)
		funct3.plot(theta[i],omega[i])

		if i%100 ==0:
				funct4.plot(theta[i],omega[i])

	if i>50000:
		rate (3000)

		if i%100 ==0:
				funct4.plot(theta[i],omega[i])

