# Chaotic-system
This is the final project I worked on for PHYS304 Computational Physics with Professor Daniel Grin at Haverford College in Spring 2020 semester.
This project studied two Chaotic systems, both of which are governed by a system of 2nd order differential equations: Damped Driven Pendulum System and Henon-Heiles systems. This work made use of odeint from scipy.integrate to find the numerical solutions. I also explored different plots for visualizations, and paid a great focus on Poincare scetions and Bifurcation diagrams. 

Here is a breakdown of the files in this repository:
* [xia_final_Chaotic_Systems.pdf](https://github.com/shufan1/Chaotic-system/blob/master/xia_final_Chaotic_Systems.pdf): a write-up for this final project
* DDP: contains all files related to sutdying DDP systems
  * xia_DDP.py: the behavior of DDP systems with different parameters
  * xia_InitialCond.py: for comparing identical systems under different initial conditions
  * xia_Bifurcate.py: for generating bifurcation diagram and finding the critical values of bifurcation points
  * xia_animate_plot.py: a simple animation that shows the time-evolution of the system
* HenonHeiles:
  * xia_HenonHeiles.py: solve and produce 3D phase plots of Henon-Heiles systems under different inital conditions.
  * xia_Poincare.py: generate the Poincare sections of the system under different initial conditions

