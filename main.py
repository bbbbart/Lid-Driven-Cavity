import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

re = 10
t = 0.25

u = np.zeros((100,100))
v = np.zeros((100,100))
p = np.zeros((100,100))

wall_v = 1.0

ghost_u = np.zeros((102,102))
ghost_v = np.zeros((102,102))
ghost_p = np.zeros((102,102))
tolerance = 0.001

x = np.linspace(-50, 50 , 100)
y = np.linspace(-50, 50, 100)

fig, ax = plt.subplots()
fig.set_size_inches(8,6)
fig.canvas.manager.set_window_title("Lid Cavity")

splot = plt.streamplot(x, y, u, v, cmap='viridis')
contour = ax.contour(x, y, np.sqrt(v*v + u*u))

def updateGrid(frame):
    global p, u, v, contour, splot
    in_tolerance = False            

    ghost_p[1:-1, 1:-1] = p
    ghost_p[0::len(ghost_p)-1, 1:-1] =  p[0::len(p)-1]
    ghost_p[1:-1, 0::len(ghost_p)-1] = p[:, 0::len(p)-1]

    ghost_v[1:-1, 1:-1] = v
    ghost_v[0::len(ghost_v)-1, 1:-1] = -v[0::len(v)-1]
    ghost_v[1:-1, 0::len(ghost_v)-1] = -v[:, 0::len(v)-1]
    
    ghost_u[1:-1, 1:-1] = u
    ghost_u[0, 1:len(ghost_u)-1] =  -u[0,:] + 2*wall_v
    ghost_u[1:-1, 0::len(ghost_u)-1] = -u[:, 0::len(u)-1]

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

    old_v = v
    old_u = u
    
    u = old_u + t*(-u*(left_u - right_u)/2 - v*(up_u - down_u)/2) + (1/re)*(left_u + right_u + up_u + down_u - 4*u) ## horiz velocity
    v = old_v + t*(-old_u*(left_v - right_v)/2 - old_v*(up_v - down_v)/2) + (1/re)*(left_v + right_v + up_v + down_v - 4*v) ## vert velocity
    b = (1/t)*(right_u - left_u + up_v - down_v)/2 ## check error

    while not in_tolerance:
        p_new = (up_p + down_p + left_p + right_p - b)/4 ## pressure solver
        if (p_new[50,50] - p[50,50]) > tolerance:
            p = p_new

            ghost_p[0::len(ghost_p)-1, 1:-1] =  p[0::len(p)-1]
            ghost_p[1:-1, 0::len(ghost_p)-1] = p[:, 0::len(p)-1]

            ghost_p[1:-1, 1:-1] = p
            up_p = ghost_p[:-2, 1:-1] 
            down_p = ghost_p[2:, 1:-1] 
            left_p = ghost_p[1:-1, :-2] 
            right_p = ghost_p[1:-1, 2:] 

            p_new = (up_p + down_p + left_p + right_p - b)/2 ## pressure solver
        else:
            in_tolerance = True
            p = p_new

    u = u - t*(right_p - left_p)/2 ## correction
    v = v - t*(up_p - down_p)/2 ## correction
    
    ax.clear()
    contour = ax.contour(x, y, np.sqrt(v*v + u*u))
    splot = plt.streamplot(x, y, u, v, cmap='viridis')

    return(contour,splot)

ani = animation.FuncAnimation(fig, updateGrid, frames = 999, interval = 100)
plt.show()