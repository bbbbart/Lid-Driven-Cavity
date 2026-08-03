import matplotlib.widgets as widget


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
