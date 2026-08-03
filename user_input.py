import widgets as wid

def buttonPresses(event):
    global pause, running
    if event.key ==  ' ':
      pause = not pause
    elif event.key == 'r':
      resetGrid()
      running = False

def updateRe(val):
   global re
   re = wid.re_slider.val

def updateTimestep(val):
   global t
   t = wid.step_slider.val

def start(event):
    global running
    running = not running

def resetGrid():
    global u, v, p, total_v, vorticity, ghost_p, ghost_u, ghost_v
    u = np.zeros((100,100))
    v = np.zeros((100,100))
    p = np.zeros((100,100))
    total_v = np.zeros((100,100))
    vorticity = np.zeros((100,100))

    ghost_u = np.zeros((102,102))
    ghost_v = np.zeros((102,102))
    ghost_p = np.zeros((102,102))