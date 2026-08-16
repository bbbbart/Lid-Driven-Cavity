import matplotlib.animation as animation
import user_input as ui
import simulation as sim

def main():
    sim.fig.canvas.mpl_connect('key_press_event', ui.buttonPresses) 
    sim.step_slider.on_changed(ui.updateTimestep)
    sim.re_slider.on_changed(ui.updateRe)
    sim.start_button.on_clicked(ui.start)

    ani = animation.FuncAnimation(sim.fig, sim.updateGrid, frames = 999, interval = 25)
    sim.plt.show()

if __name__ == "__main__":
    main()