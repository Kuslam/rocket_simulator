from dataclasses import dataclass
from globals import *
import numpy as np

@dataclass
class RocketStage:
    """
    Store key values in a rocket stage
    """
    mass: float
    mass_prop: float
    mass_struct: float
    thrust: float
    Isp: float

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
            stage1,
            stage2,
            cd,
            area,
            v_launchpad
    ):
        # Initialize given values
        self.state = init_state

        # self.mass_prop = mass_prop
        # self.mass_struct = self.state[6] - mass_prop
        # self.Isp = Isp
        
        self.stage1 = stage1
        self.stage2 = stage2

        self.rocket_vector = init_state[0:3] / np.linalg.norm(init_state[0:3])
        self.v_launchpad = v_launchpad
        self.cd = cd
        self.area = area

        # Initialize other values
        self.stage = stage1
        self.mass_sep = stage1.mass_struct + stage2.mass
        self.state_trajectory = []
        self.extra_trajectory = [] # include extra data here
        self.extra_state = np.zeros([1,9]) # acceleration (thrust, drag, net)
        self.stage_sep = False
        self.has_fuel = True
        self.has_crashed = False
        self.too_high = False

        # Set velocity of rocket to be launchpad velocity
        self.state[3:6] = v_launchpad

    def reset_state(self, init_state):
        # States
        self.state = init_state
        self.state_trajectory = []
        self.extra_trajectory = []

        # Stage sep
        self.thrust = self.stage1.thrust
        self.Isp = self.stage1.Isp
        self.stage = self.stage1

        # Conditions
        self.has_fuel = True
        self.has_crashed = False
        self.too_high = False

    def out_of_fuel_s1(self):
        if self.state[6] > self.mass_sep:
            oof = False
        else:
            oof = True
        return oof
    
    def out_of_fuel_s2(self):
        if self.state[6] > self.stage2.mass_struct:
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
    
    def get_vel_ECEF(self):
        return self.state[3:6] - self.v_launchpad

    def update_state(self, state, state_extra=None, track_state=True):
        # Update new state
        self.state = state
        v_ECEF = self.get_vel_ECEF()
        self.rocket_vector = v_ECEF / np.linalg.norm(v_ECEF)

        # Update state if we are out of fuel/crashed
        if self.crashed_into_earth():
            self.has_crashed = True
            return

        if self.out_of_fuel():
            self.has_fuel = False
            self.state[6] = self.stage1.mass_struct
            self.stage.thrust = 0

        if self.altitutde_too_high():
            self.too_high = True
            return

        if track_state:
            self.state_trajectory.append(state)
            self.extra_trajectory.append(state_extra)
        
        return