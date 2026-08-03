def main():
    fig.canvas.mpl_connect('key_press_event', buttonPresses)
    step_slider.on_changed(updateTimestep)
    re_slider.on_changed(updateRe)
    start_button.on_clicked(start)

    ani = animation.FuncAnimation(fig, updateGrid, frames = 999, interval = 25)
    plt.show()

if __name__ == "__main__":
    main()