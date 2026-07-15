import numpy as np
import matplotlib.pyplot as plt

re = 1.0
t = 1
u = np.zeros((101,101))
v = np.zeros((101,101))
p = np.zeros((101,101))
ghost_u = np.zeros((103,103))
ghost_v = np.zeros((103,103))
ghost_p = np.zeros((103,103))
tolerance = 0.001

ghost_u = u[1:-1, 1:-1]
ghost_v = v[1:-1, 1:-1]
ghost_p = p[1:-1, 1:-1]

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

u = u + t(-u((left_u - right_u))/2 - v((up_u - down_u)/2)) + (1/re)(left_u + right_u + up_u + down_u - 4*u) ## horiz velocity
v = u + t(-u((left_v - right_v))/2 - v((up_v - down_v)/2)) + (1/re)(left_v + right_v + up_v + down_v - 4*u) ## vert velocity

b = (right_u - left_u + up_v - down_v)/2 ## check error

p = (up_p + down_p + left_p + right_p - b)/4 ## pressure solver

u = u - (right_p - left_p)/2 ## correction
v = v - (up_p - down_p)/2 ## correction