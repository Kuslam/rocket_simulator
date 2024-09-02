from globals import R_EARTH

import numpy as np
import matplotlib.pyplot as plt

def create_plots(state_trajectory, t, extra_trajectory=None, filepath=None):
    rx = [state[0]/1000 for state in state_trajectory]
    ry = [state[1]/1000 for state in state_trajectory]
    rz = [state[2]/1000 for state in state_trajectory]
    vx = [state[3]/1000 for state in state_trajectory]
    vy = [state[4]/1000 for state in state_trajectory]
    vz = [state[5]/1000 for state in state_trajectory]
    m = [state[6]/1000 for state in state_trajectory]

    ax_thrust = [state[0] for state in extra_trajectory]
    ay_thrust = [state[1] for state in extra_trajectory]
    az_thrust = [state[2] for state in extra_trajectory]
    ax_drag = [state[3] for state in extra_trajectory]
    ay_drag = [state[4] for state in extra_trajectory]
    az_drag = [state[5] for state in extra_trajectory]
    ax_net = [state[6] for state in extra_trajectory]
    ay_net = [state[7] for state in extra_trajectory]
    az_net = [state[8] for state in extra_trajectory]
    ax_grav = [state[9] for state in extra_trajectory]
    ay_grav = [state[10] for state in extra_trajectory]
    az_grav = [state[11] for state in extra_trajectory]
    flight_path_angle = [state[12] for state in extra_trajectory]
    altitude = [state[13] for state in extra_trajectory]

    # Create a plot
    plt.figure(figsize=(15, 8))

    # Plot positions
    plt.subplot(2, 5, 1)
    plt.plot(t, rx, 'r-', label='rx [km]', linewidth=1.5)
    plt.plot(t, ry, 'g-', label='ry [km]', linewidth=1.5)
    plt.plot(t, rz, 'b-', label='rz [km]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Positions [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot position changes
    plt.subplot(2, 5, 2)
    plt.plot(t, rx-rx[0], 'r-', label='rx [km]', linewidth=1.5)
    plt.plot(t, ry-ry[0], 'g-', label='ry [km]', linewidth=1.5)
    plt.plot(t, rz-rz[0], 'b-', label='rz [km]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Changes in Position from Start [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot velocities
    plt.subplot(2, 5, 3)
    plt.plot(t, vx, 'r-', label='vx [km/s]', linewidth=1.5)
    plt.plot(t, vy, 'g-', label='vy [km/s]', linewidth=1.5)
    plt.plot(t, vz, 'b-', label='vz [km/s]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Velocities [ECI]')
    plt.legend()
    plt.grid(True)

    # Plot mass
    plt.subplot(2, 5, 4)
    plt.plot(t, m, 'r-', label='m [kg]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Mass')
    plt.ylabel('x 1e3')
    plt.legend()
    plt.grid(True)

    # Plot net accelerations
    plt.subplot(2, 5, 5)
    plt.plot(t, ax_net, 'r-', label='ax [km/s^2]', linewidth=1.5)
    plt.plot(t, ay_net, 'g-', label='ay [km/s^2]', linewidth=1.5)
    plt.plot(t, az_net, 'b-', label='az [km/s^2]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Net Acceleration [ECI]')
    plt.legend()
    plt.grid(True)

    # Plot thrust accelerations
    plt.subplot(2, 5, 6)
    plt.plot(t, ax_thrust, 'r-', label='ax [km/s^2]', linewidth=1.5)
    plt.plot(t, ay_thrust, 'g-', label='ay [km/s^2]', linewidth=1.5)
    plt.plot(t, az_thrust, 'b-', label='az [km/s^2]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Thrust Acceleration [ECI]')
    plt.legend()
    plt.grid(True)

    # Plot drag accelerations
    plt.subplot(2, 5, 7)
    plt.plot(t, ax_drag, 'r-', label='ax [km/s^2]', linewidth=1.5)
    plt.plot(t, ay_drag, 'g-', label='ay [km/s^2]', linewidth=1.5)
    plt.plot(t, az_drag, 'b-', label='az [km/s^2]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Drag Acceleration [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot gravity accelerations
    plt.subplot(2, 5, 8)
    plt.plot(t, ax_grav, 'r-', label='ax [km/s^2]', linewidth=1.5)
    plt.plot(t, ay_grav, 'g-', label='ay [km/s^2]', linewidth=1.5)
    plt.plot(t, az_grav, 'b-', label='az [km/s^2]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Gravity Acceleration [ECEF]')
    plt.legend()
    plt.grid(True)

    # Plot flight path angle
    plt.subplot(2, 5, 9)
    plt.plot(t, flight_path_angle, 'r-', label='$\phi$ [deg]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Flight Path Angle')
    plt.legend()
    plt.grid(True)

    # Plot altitude
    plt.subplot(2, 5, 10)
    plt.plot(t, altitude, 'r-', label='h [km]', linewidth=1.5)
    plt.xlabel('Time [s]')
    plt.title('Altitude')
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.subplots_adjust(hspace=0.4)
    if filepath != None:
        plt.savefig(f"{filepath}/trajectory_plots.png",dpi=600)
    plt.show()
    