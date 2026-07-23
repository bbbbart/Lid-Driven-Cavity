import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import njit

re = 10
t = 0.005
wall_v = 1.0
tolerance = 0.0001
max_i = 2500
delta_x = 1.0
delta_y = 1.0
h = delta_x

u = np.zeros((100,100))
v = np.zeros((100,100))
p = np.zeros((100,100))

ghost_u = np.zeros((102,102))
ghost_v = np.zeros((102,102))
ghost_p = np.zeros((102,102))

x = np.linspace(50, -50, 100)
y = np.linspace(50, -50, 100)

fig, ax = plt.subplots()
fig.set_size_inches(8,6)
fig.canvas.manager.set_window_title("Lid Cavity")

quiv = ax.quiver(x[::10], y[::10], u[::10, ::10], v[::10, ::10], scale=15)

@njit
def pressureSolver(grid, ghost, b_):
    for i in range(max_i):
        error = 0
        ghost[0, 1:-1] =  grid[0, :]
        ghost[-1, 1:-1] =  grid[-1, :]
        ghost[1:-1, -1] = grid[:, -1]
        ghost[1:-1, 0] = grid[:, 0]
        grid = ghost[1:-1, 1:-1] 
        for y in range(1, 101):
            for x in range(1, 101):
                old_p = ghost[y,x]

                ghost[y,x] = (ghost[y + 1, x] + ghost[y - 1, x] + ghost[y, x - 1] + ghost[y, x + 1] - (h**2)*b_[y-1,x-1])/4 ## pressure solver
            
                error = max(error, abs(ghost[y,x]-old_p))

        if (error < tolerance):
            return ghost, grid
            break
    return ghost, grid


def updateGrid(frame):
    for i in range(250):
        global p, u, v, quiv, ghost_p
        u[0,:] = 1.0
        v[0,:] = 0.0

        ghost_v[1:-1, 1:-1] = v
        ghost_v[[0, -1], 1:-1] = -v[[0, -1], :]
        ghost_v[1:-1, [0, -1]] = -v[:, [0, -1]]

        up_v = ghost_v[:-2, 1:-1] 
        down_v = ghost_v[2:, 1:-1] 
        left_v = ghost_v[1:-1, :-2] 
        right_v = ghost_v[1:-1, 2:] 

        ghost_u[1:-1, 1:-1] = u
        ghost_u[0, 1:-1] =  -1*(u[0,:]) + 2*wall_v
        ghost_u[-1, 1:-1] = -u[-1,:]
        ghost_u[1:-1, [0, -1]] = -u[:, [0, -1]]

        up_u = ghost_u[:-2, 1:-1] 
        down_u = ghost_u[2:, 1:-1] 
        left_u = ghost_u[1:-1, :-2] 
        right_u = ghost_u[1:-1, 2:]

        old_v = v.copy()
        old_u = u.copy()

        u = old_u + t*((-old_u*(right_u - left_u )/(2*delta_x) - old_v*(down_u - up_u)/(2*delta_y)) + (1/re)*(((left_u - 2*old_u + right_u)/delta_x**2)+((up_u - 2*old_u + down_u)/delta_y**2))) ## horiz velocity
        v = old_v + t*((-old_u*(right_v - left_v)/(2*delta_x) - old_v*(down_v - up_v)/(2*delta_y)) + (1/re)*(((left_v - 2*old_v + right_v)/delta_x**2)+((up_v - 2*old_v + down_v)/delta_y**2))) ## horiz velocity

        u[0,:] = 1.0
        v[0,:] = 0.0

        ghost_v[1:-1, 1:-1] = v
        ghost_v[[0, -1], 1:-1] = -v[[0, -1], :]
        ghost_v[1:-1, [0, -1]] = -v[:, [0, -1]]

        ghost_u[1:-1, 1:-1] = u
        ghost_u[0, 1:-1] =  -1*(u[0,:]) + 2*wall_v
        ghost_u[-1, 1:-1] = -u[-1,:]
        ghost_u[1:-1, [0, -1]] = -u[:, [0, -1]]        

        left_u = ghost_u[1:-1, :-2] 
        right_u = ghost_u[1:-1, 2:] 

        up_v = ghost_v[:-2, 1:-1] 
        down_v = ghost_v[2:, 1:-1] 

        b = (1/t)*((right_u - left_u)/(2*delta_x) + (down_v - up_v)/(2*delta_y)) ## check error

        ghost_p, p = pressureSolver(p, ghost_p, b)

        up_p = ghost_p[:-2, 1:-1] 
        down_p = ghost_p[2:, 1:-1] 
        left_p = ghost_p[1:-1, :-2] 
        right_p = ghost_p[1:-1, 2:] 

        u = u - t*(right_p - left_p)/(2*delta_x) ## correction
        v = v - t*(down_p - up_p)/(2*delta_y) ## correction

    quiv.set_UVC(u[::10, ::10], v[::10, ::10])
    return(quiv)

ani = animation.FuncAnimation(fig, updateGrid, frames = 999, interval = 25)
plt.show()