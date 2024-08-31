import numpy as np
import matplotlib.pyplot as plt

def create_plots(state_trajectories, t):
    rx = [state[0]/1000 for state in state_trajectories]
    ry = [state[1]/1000 for state in state_trajectories]
    rz = [state[2]/1000 for state in state_trajectories]
    vx = [state[3]/1000 for state in state_trajectories]
    vy = [state[4]/1000 for state in state_trajectories]
    vz = [state[5]/1000 for state in state_trajectories]
    m = [state[6] for state in state_trajectories]
    print(m)

    # Create a plot
    plt.figure(figsize=(12, 6))

    # Plot positions
    plt.subplot(3, 1, 1)
    plt.plot(t, rx, 'r-', label='rx [km]', linewidth=1.5)
    plt.plot(t, ry, 'g-', label='ry [km]', linewidth=1.5)
    plt.plot(t, rz, 'b-', label='rz [km]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Positions [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot velocities
    plt.subplot(3, 1, 2)
    plt.plot(t, vx, 'r-', label='vx [km/s]', linewidth=1.5)
    plt.plot(t, vy, 'g-', label='vy [km/s]', linewidth=1.5)
    plt.plot(t, vz, 'b-', label='vz [km/s]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocties [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot mass
    plt.subplot(3, 1, 3)
    plt.plot(t, m, 'r-', label='vx [km/s]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Mass [kg]')
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.subplots_adjust(hspace=0.4)
    plt.show()
    