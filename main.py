import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.widgets as widget
from numba import jit

re = 500
t = 0.1
wall_v = 1.0
tolerance = 0.001
max_i = 2500
delta_x = 1.0
delta_y = 1.0
h = delta_x
pause = False

u = np.zeros((100,100))
v = np.zeros((100,100))
p = np.zeros((100,100))
total_v = np.zeros((100,100))
vorticity = np.zeros((100,100))

ghost_u = np.zeros((102,102))
ghost_v = np.zeros((102,102))
ghost_p = np.zeros((102,102))

x = np.linspace(50, -50 , 100)
y = np.linspace(50, -50, 100)

fig, ax = plt.subplots(1, 3, sharex=True, sharey=True)
fig.set_size_inches(20, 7)
fig.canvas.manager.set_window_title("Velocity Quiver Plot")
ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

plt.tight_layout()

quiv = ax[0].quiver(x[::5], y[::5], u[::5, ::5], v[::5, ::5], scale=8, pivot='mid')
v_map = ax[1].imshow(np.fliplr(total_v), vmin=0, vmax=0.5, cmap='viridis', extent=[-50, 50, -50, 50], interpolation = 'bilinear')
vort_map = ax[2].imshow(np.fliplr(vorticity), vmin=-0.1, vmax=0.1, cmap='jet', extent=[-50, 50, -50, 50], interpolation = 'bilinear')

def buttonPresses(event):
    global pause
    if event.key ==  ' ':
      pause = not pause

fig.canvas.mpl_connect('key_press_event', buttonPresses)

@jit
def pressureSolver(grid, ghost, b_):
        for i in range(max_i):
            error = 0
            ghost[0, 1:-1] =  grid[0, :]
            ghost[-1, 1:-1] =  grid[-1, :]
            ghost[1:-1, -1] = grid[:, -1]
            ghost[1:-1, 0] = grid[:, 0]
            for y in range(1, 101):
                for x in range(1, 101):
                    old_p = ghost[y,x]

                    g = (ghost[y + 1, x] + ghost[y - 1, x] + ghost[y, x - 1] + ghost[y, x + 1] - (h**2)*b_[y-1,x-1])/4 ## pressure solver
                    ghost[y,x] = ghost[y,x] + 1.7*(g - ghost[y,x])

                    error = max(error, abs(ghost[y,x]-old_p))
            if (error < tolerance):
                ghost[0, 1:-1] =  grid[0, :]
                ghost[-1, 1:-1] =  grid[-1, :]
                ghost[1:-1, -1] = grid[:, -1]
                ghost[1:-1, 0] = grid[:, 0]
                return ghost, grid
                break
        ghost[0, 1:-1] = grid[0, :]
        ghost[-1, 1:-1] = grid[-1, :]
        ghost[1:-1, -1] = grid[:, -1]
        ghost[1:-1, 0] = grid[:, 0]

        return ghost, grid

def updateGrid(frame):
    if not pause:
        for i in range(100):
            global p, u, v, quiv, ghost_p, total_v, vorticity
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

            u = old_u + t*(-old_u*(np.where(old_u >= 0, (old_u - left_u )/delta_x, (right_u - old_u)/delta_x)) - old_v*(np.where(old_v >= 0, (old_u - up_u)/delta_y, (down_u - old_u)/delta_y)) + (1/re)*(((left_u - 2*old_u + right_u)/delta_x**2) + ((up_u - 2*old_u + down_u)/delta_y**2))) ## horiz velocity
            v = old_v + t*(-old_u*(np.where(old_u >= 0, (old_v - left_v )/delta_x, (right_v - old_v)/delta_x)) - old_v*(np.where(old_v >= 0, (old_v - up_v)/delta_y, (down_v - old_v)/delta_y)) + (1/re)*(((left_v - 2*old_v + right_v)/delta_x**2) + ((up_v - 2*old_v + down_v)/delta_y**2))) ## vert velocity

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

            p = ghost_p[1:-1, 1:-1]
            up_p = ghost_p[:-2, 1:-1] 
            down_p = ghost_p[2:, 1:-1] 
            left_p = ghost_p[1:-1, :-2] 
            right_p = ghost_p[1:-1, 2:]

            u = u - t*(right_p - left_p)/(2*delta_x) ## correction
            v = v - t*(down_p - up_p)/(2*delta_y) ## correction

            u[0,:] = 1.0
            v[0,:] = 0.0

            up_v = ghost_v[:-2, 1:-1] 
            down_v = ghost_v[2:, 1:-1] 
            left_u = ghost_u[1:-1, :-2] 
            right_u = ghost_u[1:-1, 2:]

            vorticity = ((up_v - down_v)/2*delta_x) - ((right_v - left_v)/2*delta_y)

        ax[0].set_title("Quiver")        
        ax[1].set_title("Total Velocity")
        ax[2].set_title("Vorticity")

        total_v = np.sqrt(u**2 + v**2)
        quiv.set_UVC(u[::5, ::5], v[::5, ::5])
        v_map.set_data(np.fliplr(total_v))
        vort_map.set_data(np.fliplr(vorticity))
        return(v_map, quiv)

    else:
        ax[0].set_title("Quiver (Paused)")        
        ax[1].set_title("Total Velocity (Paused)")
        ax[2].set_title("Vorticity (Paused)")

        quiv.set_UVC(u[::5, ::5], v[::5, ::5])
        v_map.set_data(np.fliplr(total_v))
        vort_map.set_data(np.fliplr(vorticity))
        return(v_map, quiv)
ani = animation.FuncAnimation(fig, updateGrid, frames = 999, interval = 25)
plt.show()