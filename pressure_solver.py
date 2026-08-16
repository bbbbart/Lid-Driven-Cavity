from numba import jit
from variables import max_i, grid_size, tolerance, h

@jit
def pressureSolver(grid, ghost, b_):
        for i in range(max_i):
            error = 0

            ## set boundary conditions

            ghost[0, 1:-1] =  grid[0, :]
            ghost[-1, 1:-1] =  grid[-1, :]
            ghost[1:-1, -1] = grid[:, -1]
            ghost[1:-1, 0] = grid[:, 0]

            for y in range(1, grid_size + 1):
                for x in range(1, grid_size + 1):
                    old_p = ghost[y,x]

                    ##  pressure Poisson equation

                    g = (ghost[y + 1, x] + ghost[y - 1, x] + ghost[y, x - 1] + ghost[y, x + 1] - (h**2)*b_[y-1,x-1])/4

                    ## SOR multiplier

                    ghost[y,x] = ghost[y,x] + 1.7*(g - ghost[y,x])

                    ## check error

                    error = max(error, abs(ghost[y,x]-old_p))
            if (error < tolerance):

                ## set boundary conditions

                ghost[0, 1:-1] =  grid[0, :]
                ghost[-1, 1:-1] =  grid[-1, :]
                ghost[1:-1, -1] = grid[:, -1]
                ghost[1:-1, 0] = grid[:, 0]
                return ghost, grid
                break

        ## set boundary conditions

        ghost[0, 1:-1] = grid[0, :]
        ghost[-1, 1:-1] = grid[-1, :]
        ghost[1:-1, -1] = grid[:, -1]
        ghost[1:-1, 0] = grid[:, 0]
        return ghost, grid