# Lid Driven Cavity

This project uses Matplotlib to visualize a finite difference implementation of a lid driven cavity using the Navier Stokes equations and the projection method.

## Features

- 100x100 grid-based fluid simulation
- Adjustable Reynolds number and timestep
- Pressure solver using successive over-relaxation
- Heat maps showing total velocity and vorticity
- Quiver plot showing the direction of flow

### Navier Stokes Equation

The Navier Stokes equation is used to model the motion of fluids. The most basic dimensionless form of the equation is:

### $\frac{\partial u}{\partial t}+(u\cdot \nabla )u= -\nabla  p+\frac{1}{Re}\nabla ^{2}u$

u = velocity vector
p = pressure
Re = Reynolds number
t = time

This simulation uses this equation with the projection method. Using this method, an initial velocity is calculated, the pressure Poisson equation is used to calculate the pressure field, and finally the pressure field is used to adjust the velocity.

#### Initial Velocity

At each timestep, the initial velocities of the grid can be calculated with a form of the Navier Stokes equation without pressure values:

### x: $\frac{\partial u}{\partial t}=-u\frac{\partial u}{\partial x}-v\frac{\partial u}{\partial y}+\frac{1}{Re}(\frac{\partial ^{2}u}{\partial x^{2}}+\frac{\partial ^{2}u}{\partial y^{2}})$

### y: $\frac{\partial v}{\partial t}=-u\frac{\partial v}{\partial x}-v\frac{\partial v}{\partial y}+\frac{1}{Re}(\frac{\partial ^{2}v}{\partial x^{2}}+\frac{\partial ^{2}v}{\partial y^{2}})$

u = horizontal velocity
v = vertical velocity

Equations that use partial differentials need to be put into a finite difference form that approximates derivatives. Earlier versions of the project used central difference, which takes the average of all neighboring cells. The current version uses a more stable method called first-order upwind, which uses conditional statements to check which direction the fluid is flowing. These statements then decide which of the cell's neighbors should be used in the finite difference approximation. For example a positive horizontal velocity would have the approximation $\frac{u_{i,j}-u_{i-1,j}}{\Delta x}$. A negative velocity would have the approximation of $\frac{u_{i+1,j}-u_{i,j}}{\Delta x}$.

### Pressure

In a simulation, the amount of fluid leaving must be equal to the amount entering. This can be written mathematically as:

### $\nabla \cdot u = 0$

Continuity is maintained by using the pressure Poisson equation, which creates a pressure field that is used to adjust the velocity. After calculating a temporary velocity, the divergence from continuity can be found by the equation:

### $b_{i,j}=\frac{1}{\Delta t}\nabla \cdot u^{*}$

where b is the source term which is used to calculate the pressure field.

Once b is calculated, pressure can be adjusted by repeatedly solving the finite difference equation:

### $p_{i,j} = \frac{p_{i+1,j} +p_{i-1,j}+p_{i,j+1}+p_{i,j-1}-h^{2}\cdot b_{i,j}}{4}$

This equation adjusts the pressure field to better maintain continuity. The equation repeats until the pressure difference between one pressure iteration and another is within a tolerance of 0.001.

To speed up the solver, this simulation uses a form of the Gauss-Seidel method called successive over-relaxation. In Gauss Seidel, once a cell in the grid is calculated, it is immediately used in the next calculation. The simulation does this by using a for loop that updates one cell at a time. Successive over-relaxation uses a sort of multiplier on the newly calculated pressure. This causes the pressure to jump ahead each iteration, closer to the final pressure field. The default multiplier value is 1.7.

### Velocity Correction

After solving for the pressure field, the velocity is adjusted using:

### x: $u_{i,j}^{n+1}=u_{i,j}^{n+1}-\Delta t(\frac{p_{i+1,j}-p_{i-1,j}}{2\Delta x})$
### y: $v_{i,j}^{n+1}=v_{i,j}^{n+1}-\Delta t(\frac{p_{i,j+1}-p_{i,j-1}}{2\Delta y})$

### Boundary Conditions

To maintain boundary conditions and for use in calculations, the simulation uses a ghost grid which is larger than the simulation grid. On a top moving wall, a Dirichlet (no-slip) condition sets the cells on the border equal to the wall velocity. For stationary walls, another no-slip condition sets the velocity on the wall to zero.

### Grids

The window shows 3 different grids.

#### Quiver Plot: 
A grid of arrows that shows the direction of fluid flow.

#### Total Velocity: 
A heatmap showing the magnitude of the total velocity at each cell.

#### Vorticity: 
A heatmap showing the vorticity at each cell. Vorticity measures local rotations that happen at velocity gradients (shear flow).

### UI

At launch, two sliders can be adjusted before clicking the start button.

Re slider - Can be set from 10 to 5000. Adjusts the Reynolds number, which is the ratio between inertial forces and viscous forces. Low values typically create a large central vortex, while large values typically create more complex flow with small sub-vortices.
Timestep slider = Can be set from 0.05 to 1. Timestep affects the speed of the simulation by changing how much time passes each step. This can also influence stability. Lower values are typically more stable than high values.

When a simulation is running, press `space` to pause and click `r` to reset and go back to the setup menu. If the plots suddenly disappear, it is likely because the numerical solver became unstable. Reset the grid and choose different setup values.

## Requirements and Running

Download the latest release and run the executable `Lid.Driven.Cavity.exe`.

or

Clone the repository and run `pip install -r requirements.txt` to install all required libraries.

## More Info
