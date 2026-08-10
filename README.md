# Lid Cavity

This project uses MatplotLib to visualize a finite difference implementation of a lid cavity using the Navier Stokes equations and the projection method.

## Features

- 100x100 grid-based fluid simulation
- Adjustable Reynolds number and timestep
- Pressure solver using successive over-relaxation
- Heat map showing total velocity and vorticity
- Quiver plot showing the direction of flow

### Navier Stokes Equation

The Navier Stokes equation is used to model the motion of fluids. The most basic 1D dimensionless form of the equation is:

### $\frac{\partial u}{\partial t}+(u\cdot \nabla )u= -\nabla  p+\frac{1}{Re}\nabla ^{2}u$

u = velocity vector
p = pressure
Re = reynolds number
t = time

This simulation uses this equation with the projection method. Using this method, an initial velocity is calculated, the pressure poison equation is used to calculate the pressure field, and finally the pressure field is used to adjust the velocity.

#### Initial Velocity

The initial velocity of the grid can be calculated with a form of the navier stokes equation without pressure. The finite difference form in the x and y directions are:

### x: $u_{i,j}=u_{i,j}+\Delta t(-u_{i,j}\frac{u_{i+1,j}-u_{i-1,j}}{2\Delta x}-v_{i,j}\frac{u_{i,j}-u_{i,j}}{2\Delta y})+\frac{1}{Re}(\frac{u_{i,j} - 2u_{i,j}+u_{i,j}}{x^{2}}+\frac{u_{i,j} - 2u_{i,j}+u_{i,j}}{y^{2}})$

### y: $v_{i,j}=v_{i,j}+\Delta t(-u_{i,j}\frac{v_{i+1,j}-v_{i-1,j}}{2\Delta x}-v_{i,j}\frac{v_{i,j}-v_{i,j}}{2\Delta y})+\frac{1}{Re}(\frac{v_{i,j} - 2v_{i,j}+v_{i,j}}{x^{2}}+\frac{v_{i,j} - 2v_{i,j}+v_{i,j}}{y^{2}})$

u = horizontal velocity
v = vertical velocity

### Pressure

In a simulation, the amount of fluid leaving must be equal to the amount entering. This can be written mathematically as:

### $\nabla \cdot u = 0$

Continuity is maintained by using the pressure poison equation, which adjusts velocity to control the amount of fluid crossing the simulation boundary. After calculating a temporary velocity, the divergence from continuity can be found by the equation:

### $b_{i,j}=\frac{1}{\Delta t}\Delta \cdot u^{*}$

where $b_{i,j} =\nabla ^{2}p$

Once b is calculated, pressure can be adjusted by repeatedly solving the finite difference equation:

### $p_{i,j} = \frac{p_{i+1,j} +p_{i-1,j}+p_{i,j+1}+p_{i,j-1}-h^{2}\cdot b_{i,j}}{4}$

This equation adjusts the pressure field to better maintain continuity. The equation repeats until the pressure difference between one timestep and another is within a tolerance of 0.001.

### Velocity Correction



### Boundary Conditions

### UI

## Requirements and Running

## More Info
