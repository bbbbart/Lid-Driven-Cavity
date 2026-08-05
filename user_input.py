import simulation as sim
import variables as var

def buttonPresses(event):
    if event.key ==  ' ':
      sim.pause = not sim.pause
    elif event.key == 'r':
      resetGrid()
      sim.running = False

def updateRe(val):
   var.re = sim.re_slider.val

def updateTimestep(val):
   var.t = sim.step_slider.val

def start(event):
    sim.running = not sim.running

def resetGrid():
    sim.u = sim.np.zeros((var.grid_size, var.grid_size))
    sim.v = sim.np.zeros((var.grid_size, var.grid_size))
    sim.p = sim.np.zeros((var.grid_size, var.grid_size))
    sim.total_v = sim.np.zeros((var.grid_size, var.grid_size))
    sim.vorticity = sim.np.zeros((var.grid_size, var.grid_size))

    sim.ghost_u = sim.np.zeros((var.grid_size + 2, var.grid_size + 2))
    sim.ghost_v = sim.np.zeros((var.grid_size + 2, var.grid_size + 2))
    sim.ghost_p = sim.np.zeros((var.grid_size + 2, var.grid_size + 2))