import numpy as np
import visuals as vis

# Define global variables
G0 = 9.81 #m/s
MU = 3.986e14 #m^3/s^2
R_EARTH = 6350e3 #m

class Rocket:
    """
    Rocket object contains key information about the rocket including:
    - state: [rx, ry, rz, vx, vy, vz, m]
      - rx, ry, rz: position in ECEF (in m)
      - vx, vy, vz: velocity in the x, y, and z directions (in m/s)
      - m: mass (in kg)
    - mass_prop: mass of propellant (in kg)
    - Isp: Isp of rocket (in s)
    - thrust: thrust of rocket (in N)
    - rocket_vector: unit vector of direction rokcet is pointing (in ECI)
    """
    def __init__(
            self,
            init_state,
            mass_prop,
            Isp,
            thrust,
            cd,
            area,
            v_launchpad
    ):
        # Initialize given values
        self.state = init_state
        self.mass_prop = mass_prop
        self.mass_struct = self.state[6] - mass_prop
        self.Isp = Isp
        self.thrust = thrust
        self.rocket_vector = init_state[0:3] / np.linalg.norm(init_state[0:3])
        self.v_launchpad = v_launchpad
        self.cd = cd
        self.area = area

        # Initialize other values
        self.state_trajectory = []
        self.extra_trajectory = [] # include extra data here
        self.extra_state = np.zeros([1,9]) # acceleration (thrust, drag, net)
        self.has_fuel = True
        self.has_crashed = False
        self.too_high = False

        # Set velocity of rocket to be launchpad velocity
        self.state[3:6] = v_launchpad

    def reset_state(self, init_state):
        self.state = init_state
        self.state_trajectory = []
        self.extra_trajectory = []
        self.has_fuel = True
        self.has_crashed = False
        self.too_high = False

    def out_of_fuel(self):
        if self.state[6] > self.mass_struct:
            oof = False
        else:
            oof = True
        return oof
    
    def crashed_into_earth(self):
        if self.get_altitude() <= R_EARTH:
            crashed = True
        else:
            crashed = False
        return crashed

    def altitutde_too_high(self):
        if self.get_altitude() >= 368000e3:
            too_high = True
        else:
            too_high = False
        return too_high

    def get_altitude(self):
        return np.linalg.norm(self.state[0:3])

    def update_state(self, state, state_extra=None, track_state=True):
        # Update new state
        self.state = state
        v_ECEF = state[3:6] - self.v_launchpad
        self.rocket_vector = v_ECEF / np.linalg.norm(v_ECEF)

        # Update state if we are out of fuel/crashed
        if self.crashed_into_earth():
            self.has_crashed = True

        if self.out_of_fuel():
            self.has_fuel = False
            self.state[6] = self.mass_struct
            self.thrust = 0

        if self.altitutde_too_high():
            self.too_high = True

        if track_state:
            self.state_trajectory.append(state)
            self.extra_trajectory.append(state_extra)

def create_rocket():
    # Initialize state
    # rx, ry, rz = 916.361e3, -5539.907e3, 3014.812e3 #m
    rx, ry, rz = 0, 0, 6350e3
    vx, vy, vz = 0, 0, 0 #m/s
    mass = 549000 #kg
    mass_struct = 25600 #kg
    mass_prop = mass - mass_struct #kg
    Isp = 292 #s
    thrust = 7600e3 #N
    cd = 0.75
    area = np.pi/4 * 3.7**2

    # Assume we launch out of Cape Canaveral
    v_launchpad = np.array([0, 0, 0])

    init_state = np.array([rx, ry, rz, vx, vy, vz, mass])

    rocket = Rocket(init_state=init_state,
                    mass_prop=mass_prop,
                    Isp=Isp,
                    thrust=thrust,
                    cd=cd,
                    area=area,
                    v_launchpad=v_launchpad)
    
    return rocket

def atmospheric_model(position):
    """
    Model of earth's atmosphere
    source: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    """

    altitude = np.linalg.norm(position) - R_EARTH
    if altitude < 0: altitude = 0

    # Check if over von-karman line
    if altitude > 100e3: return 0

    # Apply model
    if altitude < 11000:
        T = 15.04 - 0.00649*altitude
        p = 101.29 * ((T + 273.1)/288.08)**5.256
    elif altitude < 25000:
        T = -56.46
        p = 22.65*np.exp(1.73-0.000157*altitude)
    elif altitude:
        T = -131.21 + 0.00299*altitude
        p = 2.488 * ((T + 273.1)/216.6)**(-11.388)
    
    density = p / (0.289*(T+273.1))
    return density

def forward_step(rocket, state, dt, output_extra_state = True):
    # Get key parameters from state
    r = state[0:3]
    v = state[3:6]
    mass = state[6]

    # Find thrust
    thrust = rocket.thrust * rocket.rocket_vector

    # Find drag (opposite of rocket direction)
    v_drag = np.dot(v,rocket.rocket_vector)
    rho = atmospheric_model(r)
    drag = -0.5 * rho * rocket.area * rocket.cd * v_drag**2
    drag *= rocket.rocket_vector
    drag = np.zeros(3)
    
    # Find forces on rocket (and net acceleration)
    accel_gravity =  - MU / (np.linalg.norm(r)**3) * r
    accel_thrust = thrust / mass
    accel_drag = drag / mass
    accel_net = accel_thrust + accel_gravity + accel_drag

    # Find change in state
    state_dot = np.zeros(state.shape)
    state_dot[0:3] = v # r_dot = vel
    state_dot[3:6] = accel_net # v_dot = accel
    state_dot[6] = - np.linalg.norm(thrust) / (G0 * rocket.Isp) # mass_dot = . . . 

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

        if rocket.has_crashed:
            return t_iter
        
        if rocket.too_high:
            return t_iter

    return t_iter

def main():
    rocket = create_rocket()
    dt = 1
    t_final = 10000
    t_final = t_final - t_final%dt
    t_final_sim = simulate_trajectory(rocket=rocket,
                                      dt=dt,
                                      t_final=t_final)
    
    vis.create_plots(state_trajectory=rocket.state_trajectory,
                     t=np.linspace(0, t_final_sim, int(t_final_sim/dt)),
                     extra_trajectory=rocket.extra_trajectory)

    print("end")
    
main()