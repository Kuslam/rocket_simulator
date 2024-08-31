import numpy as np
import matplotlib.pyplot as plt

def create_plots(state_trajectories, t):
    rx = [state[0] for state in state_trajectories]
    ry = [state[1] for state in state_trajectories]
    rz = [state[2] for state in state_trajectories]
    vx = [state[3] for state in state_trajectories]
    vy = [state[4] for state in state_trajectories]
    vx = [state[5] for state in state_trajectories]
    m = [state[6] for state in state_trajectories]

    # Create a plot
    plt.figure(figsize=(10, 6))

    # Plot rx, ry, and rz over time
    plt.plot(t, rx, marker='o', linestyle='-', color='r', label='rx')
    plt.plot(t, ry, marker='s', linestyle='--', color='g', label='ry')
    plt.plot(t, rz, marker='^', linestyle='-.', color='b', label='rz')

    # Add labels and title
    plt.xlabel('Time (t)')
    plt.ylabel('Values')
    plt.title('Plot of rx, ry, and rz versus Time')
    plt.legend()

    # Show grid
    plt.grid(True)

    # Display the plot
    plt.show()
    