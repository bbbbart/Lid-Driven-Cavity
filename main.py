import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

re = 10
t = 0.01
wall_v = 1.0
tolerance = 0.01
max_i = 200
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
##splot = plt.streamplot(x, y, u, v, cmap='viridis')
##contour = ax.contour(x, y, np.sqrt(v*v + u*u))

def updateGrid(frame):
    global p, u, v, quiv

    for i in range(10):
        ghost_v[1:-1, 1:-1] = v
        ghost_v[[0, -1], 1:-1] = -v[[0, -1], :]
        ghost_v[1:-1, [0, -1]] = -v[:, [0, -1]]
        
        ghost_u[1:-1, 1:-1] = u
        ghost_u[0, 1:-1] =  -1*(u[0,:]) + 2*wall_v
        ghost_u[-1, 1:-1] = -u[-1,:]
        ghost_u[1:-1, [0, -1]] = -u[:, [0, -1]]

        up_u = ghost_u[:-2, 1:-1] 
        down_u = ghost_u[2:, 1:-1] 
        left_u = ghost_u[1:-1, :-2] 
        right_u = ghost_u[1:-1, 2:] 

        up_v = ghost_v[:-2, 1:-1] 
        down_v = ghost_v[2:, 1:-1] 
        left_v = ghost_v[1:-1, :-2] 
        right_v = ghost_v[1:-1, 2:] 

        old_v = v.copy()
        old_u = u.copy()
        
        u = old_u + t*((-old_u*(right_u - left_u )/(2*delta_x) - old_v*(down_u - up_u)/(2*delta_y)) + (1/re)*(((left_u - 2*old_u + right_u)/delta_x**2)+((up_u - 2*old_u + down_u)/delta_y**2))) ## horiz velocity
        v = old_v + t*((-old_u*(right_v - left_v)/(2*delta_x) - old_v*(down_v - up_v)/(2*delta_y)) + (1/re)*(((left_v - 2*old_v + right_v)/delta_x**2)+((up_v - 2*old_v + down_v)/delta_y**2))) ## horiz velocity
        
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
        for i in range(max_i):
            ghost_p[[0,-1], 1:-1] =  p[[0, -1], :]
            ghost_p[1:-1, [0, -1]] = p[:, [0, -1]]
            ghost_p[1:-1, 1:-1] = p

            up_p = ghost_p[:-2, 1:-1] 
            down_p = ghost_p[2:, 1:-1] 
            left_p = ghost_p[1:-1, :-2] 
            right_p = ghost_p[1:-1, 2:] 

            p_new = (up_p + down_p + left_p + right_p - (h**2)*b)/4 ## pressure solver
            if (np.max(np.abs(p_new - p)) < tolerance):
                p = p_new
                break
        p = p_new

        ghost_p[[0,-1], 1:-1] =  p[[0, -1], :]
        ghost_p[1:-1, [0, -1]] = p[:, [0, -1]]
        ghost_p[1:-1, 1:-1] = p

        up_p = ghost_p[:-2, 1:-1] 
        down_p = ghost_p[2:, 1:-1] 
        left_p = ghost_p[1:-1, :-2] 
        right_p = ghost_p[1:-1, 2:] 

        u = u - t*(right_p - left_p)/(2*delta_x) ## correction
        v = v - t*(down_p - up_p)/(2*delta_y) ## correction
    
    quiv.set_UVC(u[::10, ::10], v[::10, ::10])
    ##ax.clear()
    ##contour = ax.contour(x, y, np.sqrt(v**2 + u**2))
    ##splot = plt.streamplot(x, y, u, v, cmap='viridis')

    ##return(contour,splot)
    return(quiv)

ani = animation.FuncAnimation(fig, updateGrid, frames = 999, interval = 50)
plt.show()