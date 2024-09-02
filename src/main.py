# Custom libraries
from globals import *
from launch_dynamics import simulate_trajectory
from init import create_rocket
from rocket import Rocket, PitchOver
import visuals as vis
import filepaths as fp

# General libraries
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

def simulate_pitchover(x_pitchover):
    pitchover = PitchOver(x_pitchover=x_pitchover)
    
    rocket = create_rocket(pitchover=pitchover)
    t_final_sim = simulate_trajectory(rocket=rocket,
                                      dt=DT,
                                      t_final=T_FINAL)

    return t_final_sim, rocket

def optimization_function(x_pitchover):
    t_final_sim, rocket = simulate_pitchover(x_pitchover=x_pitchover)

    # We want to try and hit an orbit (not too specific rn)
    # To do so, the simulation should complete and we should have the lowest possible radius
    # We know two things: simulation should end & flight path angle should be 90 at the end
    cost_function = (t_final_sim-T_FINAL)**2 + (rocket.get_flight_path_angle() - 90)**2
    return cost_function

def main():
    x_bounds = opt.Bounds(lb=[0, 0, 0, 0],
                          ub=[15, 1e-8, 60, 10])
    results = opt.dual_annealing(func=optimization_function,
                                 bounds=x_bounds,
                                 maxiter=10)

    # Simulate best rocket
    x_best = results.x
    t_final_sim, rocket = simulate_pitchover(x_pitchover=x_best)

    # Plot and save data
    filepath = fp.make_data_folder()
    vis.create_plots(state_trajectory=rocket.state_trajectory,
                     t=np.linspace(0, t_final_sim, int(t_final_sim/DT)),
                     extra_trajectory=rocket.extra_trajectory,
                     filepath=filepath)
    np.savetxt(f"{filepath}/x_best.csv", x_best, delimiter=',', fmt='%s')

    return
    
main()