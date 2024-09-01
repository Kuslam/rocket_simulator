# Custom libraries
from globals import *
from launch_dynamics import simulate_trajectory
from init import create_rocket
from rocket import Rocket
import visuals as vis

# General libraries
import numpy as np

def main():
    rocket = create_rocket()
    dt = 1
    t_final = 10000
    t_final -= t_final%dt
    t_final_sim = simulate_trajectory(rocket=rocket,
                                      dt=dt,
                                      t_final=t_final)
    
    vis.create_plots(state_trajectory=rocket.state_trajectory,
                     t=np.linspace(0, t_final_sim, int(t_final_sim/dt)),
                     extra_trajectory=rocket.extra_trajectory)

    print("end")
    
main()