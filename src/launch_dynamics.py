"""
Main dynamics used during launch
"""
# Custom libraries
from globals import *
from init import create_rocket
from rocket import Rocket
import visuals as vis

# General libraries
import numpy as np

def forward_step(rocket, state, dt, output_extra_state = True):
    # Get key parameters from state
    r = state[0:3]
    v = state[3:6]
    mass = state[6]

    # Find thrust
    thrust = rocket.stage.thrust * rocket.rocket_vector

    # Find drag (opposite of rocket direction)
    v_drag = rocket.get_vel_ECEF()
    rho = atmospheric_model(r)
    drag = -0.5 * rho * rocket.area * rocket.cd * v_drag**2
    drag *= rocket.rocket_vector
    
    # Find forces on rocket (and net acceleration)
    accel_gravity =  - MU / (np.linalg.norm(r)**3) * r
    accel_thrust = thrust / mass
    accel_drag = drag / mass
    accel_net = accel_thrust + accel_gravity + accel_drag

    # Find change in state
    state_dot = np.zeros(state.shape)
    state_dot[0:3] = v # r_dot = vel
    state_dot[3:6] = accel_net # v_dot = accel
    state_dot[6] = - np.linalg.norm(thrust) / (G0 * rocket.stage.Isp) # mass_dot = . . . 

    # Get extra state data
    state_extra = np.concatenate([accel_thrust, accel_drag, accel_net, accel_gravity])

    return state_dot, state_extra

def rk4(step, state, rocket, dt):
    k1, k1e = step(rocket=rocket, state=state, dt=dt)
    k2, k2e = step(rocket=rocket, state=(state+k1*dt/2), dt = dt/2)
    k3, k3e = step(rocket=rocket, state=(state+k2*dt/2), dt = dt/2)
    k4, k4e = step(rocket=rocket, state=(state+k3*dt), dt = dt)

    state_dot = k1/6 + k2/3 + k3/3 + k4/6
    state_extra = k1e/6 + k2e/3 + k3e/3 + k4e/6
    return state_dot, state_extra

def simulate_trajectory(rocket, dt, t_final):
    # Initialize loop
    final_iter = int(np.ceil(t_final / dt))
    t_iter = 0

    # Iterate
    for i in range(final_iter):
        # Find new state
        state_dot, new_state_extra = rk4(step=forward_step,
                                         state=rocket.state,
                                         rocket=rocket,
                                         dt=dt)
        new_state = rocket.state + state_dot*dt

        # Update state
        t_iter += dt
        rocket.update_state(state=new_state,
                            state_extra=new_state_extra)

        if rocket.has_crashed or rocket.too_high:
            return t_iter-dt

    return t_iter