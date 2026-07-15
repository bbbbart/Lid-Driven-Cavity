import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

re = 1.0
t = 1

u = np.zeros((101,101))
v = np.zeros((101,101))
p = np.zeros((101,101))

ghost_u = np.zeros((103,103))
ghost_v = np.zeros((103,103))
ghost_p = np.zeros((103,103))
tolerance = 0.001
in_tolerance = False
total_velocity = np.sqrt(v*v + u*u)

x = np.linspace(-5.0, 5.0, 101)
y = np.linspace(-5.0, 5.0, 101)

fig, ax = plt.subplots()
fig.set_size_inches(8,6)
fig.canvas.manager.set_window_title("Lid Cavity")

fig = plt.streamplot(x, y, u, v, cmap='viridis')
fig = ax.contour(x, y, total_velocity)

def updateGrid(frame):
    global p

    ghost_u[1:-1, 1:-1] = u
    ghost_v[1:-1, 1:-1] = v
    ghost_p[1:-1, 1:-1] = p

    up_u = ghost_u[:-2, 1:-1] 
    down_u = ghost_u[2:, 1:-1] 
    left_u = ghost_u[1:-1, :-2] 
    right_u = ghost_u[1:-1, 2:] 

    up_v = ghost_v[:-2, 1:-1] 
    down_v = ghost_v[2:, 1:-1] 
    left_v = ghost_v[1:-1, :-2] 
    right_v = ghost_v[1:-1, 2:] 

    up_p = ghost_p[:-2, 1:-1] 
    down_p = ghost_p[2:, 1:-1] 
    left_p = ghost_p[1:-1, :-2] 
    right_p = ghost_p[1:-1, 2:] 

    u = u + t*(-u*(left_u - right_u)/2 - v*(up_u - down_u)/2) + (1/re)*(left_u + right_u + up_u + down_u - 4*u) ## horiz velocity
    v = u + t*(-u*(left_v - right_v)/2 - v*(up_v - down_v)/2) + (1/re)*(left_v + right_v + up_v + down_v - 4*u) ## vert velocity
    b = (1/t)*(right_u - left_u + up_v - down_v)/2 ## check error

    while not in_tolerance:
        if (p_new - p) > tolerance:
            p = p_new
            p_new = (up_p + down_p + left_p + right_p - b)/4 ## pressure solver
        else:
            in_tolerance = True
    u = u - t*(right_p - left_p)/2 ## correction
    v = v - t*(up_p - down_p)/2 ## correction
    
    fig.remove()
    fig = ax.contour(x, y, total_velocity)

    return(fig,)

plt.show()