# Lid Cavity

This project uses MatplotLib to visualize a finite difference implementation of a lid cavity using the Navier Stokes equations

## Features

- 100x100 grid-based fluid simulation
- Adjustable Reynolds number and timestep
- Pressure solver using successive over-relaxation
- Heat map showing total velocity and vorticity
- Quiver plot showing the direction of flow

### Navier Stokes Equations

The Navier Stokes equations are used to model the motion of fluids.

#### Continuity

In a simulation, the amount of fluid leaving must be equal to the amount entering. This can be written mathematically as:

### $\nabla \cdot u = 0$

### Pressure

Continuity is maintained by using the pressure poison equation, which adjusts velocity to control the amount of fluid crossing the simulation boundary. After using the momentum equation to calculate a temporary velocity, the divergence from continuity can be found by the equation:

### $b_{i,j}=\frac{1}{\Delta t}\Delta \cdot u^{*}$

where $b_{i,j} =\nabla ^{2}p$

Once b is calculated, pressure can be adjusted by repeatedly solving the finite difference equation:

### $p_{i,j} = \frac{p_{i+1,j} +p_{i-1,j}+p_{i,j+1}+p_{i,j-1}-h^{2}\cdot b_{i,j}}{4}$

This equation adjusts the pressure to better maintain continuity. The equation repeats until the pressure difference between one timestep and another is within a tolerance of 0.001.

### Boundary Conditions

### UI

## Requirements and Running

## More Info
