import matplotlib.widgets as widget
import matplotlib.pyplot as plt
import variables as var
import numpy as np
from pressure_solver import pressureSolver

pause = False
running = False

u = np.zeros((var.grid_size, var.grid_size))
v = np.zeros((var.grid_size, var.grid_size))
p = np.zeros((var.grid_size, var.grid_size))
total_v = np.zeros((var.grid_size, var.grid_size))
vorticity = np.zeros((var.grid_size, var.grid_size))

ghost_u = np.zeros((var.grid_size + 2, var.grid_size + 2))
ghost_v = np.zeros((var.grid_size + 2,var.grid_size + 2))
ghost_p = np.zeros((var.grid_size + 2, var.grid_size + 2))

x = np.linspace((var.grid_size/2), (-1)*(var.grid_size/2), var.grid_size)
y = np.linspace((var.grid_size/2), (-1)*(var.grid_size/2), var.grid_size)

fig, ax = plt.subplots(1, 3, sharex=True, sharey=True)
fig.set_size_inches(20, 7)
ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

ax[0].set_title("Quiver")        
ax[1].set_title("Total Velocity")
ax[2].set_title("Vorticity")

re_slider_axes = plt.axes([0.25, 0.6, 0.50, 0.08])
re_slider = widget.Slider(ax = re_slider_axes, 
                        label= 'Reynolds Number (Re)', 
                        valmin = 10, 
                        valmax = 500)

step_slider_axes = plt.axes([0.25, 0.4, 0.50, 0.08])
step_slider = widget.Slider(ax = step_slider_axes, 
                            label= 'Simulation Timestep', 
                            valmin = 0.05, 
                            valmax = 1)

start_button_axes = plt.axes([0.44, 0.18, 0.12, 0.06])
start_button = widget.Button(ax = start_button_axes, 
                        label = 'Start', 
                        color = 'lightgray', 
                        hovercolor = 'gray')

quiv = ax[0].quiver(x[::4], y[::4], u[::4, ::4], v[::4, ::4], scale=8, pivot='mid')
v_map = ax[1].imshow(np.fliplr(total_v), vmin=0, vmax=0.5, cmap='viridis', extent=[(-1)*(var.grid_size/2), var.grid_size/2, (-1)*(var.grid_size/2), var.grid_size/2], interpolation = 'bicubic')
vort_map = ax[2].imshow(np.fliplr(vorticity), vmin=-0.1, vmax=0.1, cmap='jet', extent=[(-1)*(var.grid_size/2), var.grid_size/2, (-1)*(var.grid_size/2), var.grid_size/2], interpolation = 'bicubic')

def updateGrid(frame):
    global u, v, p, total_v, vorticity, ghost_p
    if not pause and running:
        re_slider_axes.set_visible(False)
        step_slider_axes.set_visible(False)
        start_button_axes.set_visible(False)

        ax[0].set_visible(True)        
        ax[1].set_visible(True)        
        ax[2].set_visible(True) 
        for i in range(50):
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
            ghost_u[0, 1:-1] =  -1*(u[0,:]) + 2*var.wall_v
            ghost_u[-1, 1:-1] = -u[-1,:]
            ghost_u[1:-1, [0, -1]] = -u[:, [0, -1]]

            up_u = ghost_u[:-2, 1:-1] 
            down_u = ghost_u[2:, 1:-1] 
            left_u = ghost_u[1:-1, :-2] 
            right_u = ghost_u[1:-1, 2:]

            old_v = v.copy()
            old_u = u.copy()

            u = old_u + var.t*(-old_u*(np.where(old_u >= 0, (old_u - left_u )/var.delta_x, (right_u - old_u)/var.delta_x)) - old_v*(np.where(old_v >= 0, (old_u - up_u)/var.delta_y, (down_u - old_u)/var.delta_y)) + (1/var.re)*(((left_u - 2*old_u + right_u)/var.delta_x**2) + ((up_u - 2*old_u + down_u)/var.delta_y**2))) ## horiz velocity
            v = old_v + var.t*(-old_u*(np.where(old_u >= 0, (old_v - left_v )/var.delta_x, (right_v - old_v)/var.delta_x)) - old_v*(np.where(old_v >= 0, (old_v - up_v)/var.delta_y, (down_v - old_v)/var.delta_y)) + (1/var.re)*(((left_v - 2*old_v + right_v)/var.delta_x**2) + ((up_v - 2*old_v + down_v)/var.delta_y**2))) ## vert velocity

            u[0,:] = 1.0
            v[0,:] = 0.0

            ghost_v[1:-1, 1:-1] = v
            ghost_v[[0, -1], 1:-1] = -v[[0, -1], :]
            ghost_v[1:-1, [0, -1]] = -v[:, [0, -1]]

            ghost_u[1:-1, 1:-1] = u
            ghost_u[0, 1:-1] =  -1*(u[0,:]) + 2*var.wall_v
            ghost_u[-1, 1:-1] = -u[-1,:]
            ghost_u[1:-1, [0, -1]] = -u[:, [0, -1]]        

            left_u = ghost_u[1:-1, :-2] 
            right_u = ghost_u[1:-1, 2:] 

            up_v = ghost_v[:-2, 1:-1] 
            down_v = ghost_v[2:, 1:-1] 

            b = (1/var.t)*((right_u - left_u)/(2*var.delta_x) + (down_v - up_v)/(2*var.delta_y)) ## check error

            ghost_p, p = pressureSolver(p, ghost_p, b)

            p = ghost_p[1:-1, 1:-1]
            up_p = ghost_p[:-2, 1:-1] 
            down_p = ghost_p[2:, 1:-1] 
            left_p = ghost_p[1:-1, :-2] 
            right_p = ghost_p[1:-1, 2:]

            u = u - var.t*(right_p - left_p)/(2*var.delta_x) ## correction
            v = v - var.t*(down_p - up_p)/(2*var.delta_y) ## correction

            u[0,:] = 1.0
            v[0,:] = 0.0

            up_v = ghost_v[:-2, 1:-1] 
            down_v = ghost_v[2:, 1:-1] 
            left_u = ghost_u[1:-1, :-2] 
            right_u = ghost_u[1:-1, 2:]

            vorticity = ((up_v - down_v)/(2*var.delta_x)) - ((right_v - left_v)/(2*var.delta_y))

        ax[0].set_title("Quiver")        
        ax[1].set_title("Total Velocity")
        ax[2].set_title("Vorticity")

        total_v = np.sqrt(u**2 + v**2)
        quiv.set_UVC(u[::4, ::4], v[::4, ::4])
        v_map.set_data(np.fliplr(total_v))
        vort_map.set_data(np.fliplr(vorticity))
        return(v_map, quiv)

    elif running and pause:
        re_slider_axes.set_visible(False)
        step_slider_axes.set_visible(False)
        start_button_axes.set_visible(False)

        ax[0].set_visible(True)        
        ax[1].set_visible(True)        
        ax[2].set_visible(True) 

        ax[0].set_title("Quiver (Paused)")        
        ax[1].set_title("Total Velocity (Paused)")
        ax[2].set_title("Vorticity (Paused)")

        quiv.set_UVC(u[::4, ::4], v[::4, ::4])
        v_map.set_data(np.fliplr(total_v))
        vort_map.set_data(np.fliplr(vorticity))
        return(v_map, quiv)
    
    elif not running:
        re_slider_axes.set_visible(True)
        step_slider_axes.set_visible(True)
        start_button_axes.set_visible(True)

        ax[0].set_visible(False)        
        ax[1].set_visible(False)        
        ax[2].set_visible(False)   
    
