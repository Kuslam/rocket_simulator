# Custom libraries
from globals import *
from launch_dynamics import simulate_trajectory
from init import create_rocket
from rocket import Rocket, PitchOver
import visuals as vis

# General libraries
import numpy as np
from scipy.optimize import minimize

def simulate_pitchover(x_pitchover):
    pitchover = PitchOver(x_pitchover=x_pitchover)
    
    rocket = create_rocket(pitchover=pitchover)
    t_final_sim = simulate_trajectory(rocket=rocket,
                                      dt=DT,
                                      t_final=T_FINAL)

    return t_final_sim, rocket

def optimization_function(x_pitchover):
    t_final_sim, rocket = simulate_pitchover(x_pitchover=x_pitchover)
    return -t_final_sim

def main():
    
    # results = minimize(optimization_function, )

    # Simulate best rocket
    x_best = X0_PITCHOVER[:]
    t_final_sim, rocket = simulate_pitchover(x_pitchover=x_best)

    vis.create_plots(state_trajectory=rocket.state_trajectory,
                     t=np.linspace(0, t_final_sim, int(t_final_sim/DT)),
                     extra_trajectory=rocket.extra_trajectory)

    print("end")
    
main()