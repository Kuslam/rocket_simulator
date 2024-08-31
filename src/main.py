import numpy as np
import visuals as vis

# Define global variables
G0 = 9.81 #m/s
MU = 3.986e14 #m^3/s^2
R_EARTH = 6350e3 #m

class Rocket:
    """
    Rocket object contains key information about the rocket including:
    - state: [rx, ry, rz, vx, vy, vz, m, q0, q1, q2, q3]
      - rx, ry, rz: position in ECEF (in m)
      - vx, vy, vz: velocity in the x, y, and z directions (in m/s)
      - m: mass (in kg)
      - q0, q1, q2, q3: orientation in quaternions, with q0 being the scalar
    - mass_prop: mass of propellant (in kg)
    - Isp: Isp of rocket (in s)
    - thrust: thrust of rocket (in N)
    - thrust_vector: unit vector of direction thrust occurs in (in ECI)
    """
    def __init__(
            self,
            init_state,
            mass_prop,
            Isp,
            thrust
    ):
        self.state = init_state
        self.mass_prop = mass_prop
        self.Isp = Isp
        self.thrust = thrust
        self.thrust_vector = init_state[0:3] / np.linalg.norm(init_state[0:3]) #TODO: replace with quaternions
        self.state_trajectory = []

    def reset_state(self, init_state):
        self.state = init_state
        self.state_trajectory = []


def create_rocket():
    # Initialize state
    rx, ry, rz = 916.361e3, -5539.907e3, 3014.812e3 #m
    vx, vy, vz = 0, 0, 0 #m/s
    mass = 549000 #kg
    mass_struct = 25600 #kg
    mass_prop = mass - mass_struct #kg
    Isp = 292 #s
    thrust = 7600e3 #N

    init_state = np.array([rx, ry, rz, vx, vy, vz, mass])

    rocket = Rocket(init_state=init_state,
                    mass_prop=mass_prop,
                    Isp=Isp,
                    thrust=thrust)
    
    return rocket

def forward_step(rocket, state, dt):
    # Get key parameters from state
    r = state[0:3]
    mass = state[6]

    # See if we have thrust
    if mass > 0:
        thrust = rocket.thrust * rocket.thrust_vector
    else:
        thrust = 0 * rocket.thrust_vector
    
    # Find forces on rocket (and net acceleration)
    accel_gravity =  - MU / (np.linalg.norm(r)**3) * r
    accel_net = rocket.thrust / mass + accel_gravity

    # Find change in state
    state_dot = np.zeros(state.shape)
    state_dot[0:3] = state[3:6] # r_dot = vel
    state_dot[3:6] = accel_net # v_dot = accel
    state_dot[6] = np.linalg.norm(thrust) / (G0 * rocket.Isp) # mass_dot = . . . 

    return state_dot

def rk4(step, state, rocket, dt):
    k1 = step(rocket=rocket, state=state, dt=dt)
    k2 = step(rocket=rocket, state=(state+k1*dt/2), dt = dt/2)
    k3 = step(rocket=rocket, state=(state+k2*dt/2), dt = dt/2)
    k4 = step(rocket=rocket, state=(state+k3*dt), dt = dt)

    state_dot = k1/6 + k2/3 + k3/3 + k4/6
    return state_dot

def simulate_trajectory(rocket, dt, t_final):
    # Initialize loop
    final_iter = int(np.ceil(t_final / dt))

    # Iterate
    for i in range(final_iter):
        state_dot = rk4(step=forward_step,
                        state=rocket.state,
                        rocket=rocket,
                        dt=dt)
        rocket.state = state_dot*dt

        # Check if we ran out of fuel
        if rocket.state[6] < 0:
            rocket.state[6] = 0
        
        rocket.state_trajectory.append(rocket.state)

def main():
    rocket = create_rocket()
    dt = 0.1
    t_final = 100
    t_final = t_final - t_final%dt
    simulate_trajectory(rocket=rocket,
                        dt=dt,
                        t_final=t_final)
    
    vis.create_plots(rocket.state_trajectory, np.linspace(0, t_final, int(t_final/dt)))
    
main()