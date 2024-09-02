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

class PitchOver:
    """
    Store key values about the pitchover maneuver
    """
    def __init__(self,
                 x_pitchover):
            self.angle_inclination = np.deg2rad(x_pitchover[0])
            self.angle_azimuth = np.deg2rad(x_pitchover[1])
            self.t_start = x_pitchover[2]
            t_burn = x_pitchover[3]
            self.t_end = t_start + t_burn
        
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
            pitchover,
            cd,
            area,
            v_launchpad
    ):
        # Initialize given values
        self.state = init_state
        
        self.stage1 = stage1
        self.stage2 = stage2
        self.pitchover = pitchover

        self.rocket_vector = init_state[0:3] / np.linalg.norm(init_state[0:3])
        self.v_launchpad = v_launchpad
        self.cd = cd
        self.area = area

        # Initialize other values
        self.stage = stage1
        self.mass_sep = stage1.mass_struct + stage2.mass
        self.state_trajectory = []
        self.extra_trajectory = [] # include extra data here (e.g acceleration, flight path angle, altitude)
        self.has_separated = False
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
        self.has_crashed = False
        self.too_high = False

    def stage_separation(self):
        # To conduct stage separation
        self.stage = self.stage2
        self.state[6] = self.stage.mass
        self.has_separated = True

    def out_of_fuel(self):
        # Stage 1
        if self.state[6] < self.mass_sep and (not self.has_separated):
            oof = 1
        # Stage 2
        elif self.state[6] < self.stage2.mass_struct and self.has_separated:
            oof = 2
        else:
            oof = 0
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
    
    def get_thrust_vector(self, t):
        if t < self.pitchover.t_start or t > self.pitchover.t_end:
            # we are not doing a pitchover
            return self.rocket_vector            
        else:
            # we are doing a pitchover
            rocket_vector_rockref = np.array([0, 0, 1]) #rocket reference frame
            thrust_vector_rocketref = np.array([np.cos(self.pitchover.angle_azimuth) * np.sin(self.pitchover.angle_inclination),
                                                np.sin(self.pitchover.angle_azimuth) * np.sin(self.pitchover.angle_inclination),
                                                np.cos(self.pitchover.angle_inclination)]) #rocket reference frame
            R_rockref2ECI = find_dcm(rocket_vector_rockref, self.rocket_vector)
            thrust_vector = R_rockref2ECI @ thrust_vector_rocketref

        return thrust_vector
    
    def get_flight_path_angle(self):
        # Find the flight path angle
        r = normalize(self.state[0:3])
        v = normalize(self.state[3:6])

        # FPA is vero for IC
        if np.linalg.norm(v) == 0:
            return 0
        
        theta = np.rad2deg(np.arccos(np.dot(r,v)))
        return theta

    def update_state(self, state, state_extra=None, track_state=True):
        # Update new state
        self.state = state
        v_ECEF = self.get_vel_ECEF()
        self.rocket_vector = v_ECEF / np.linalg.norm(v_ECEF)

        # Update state if we are out of fuel/crashed
        if self.crashed_into_earth():
            self.has_crashed = True
            return

        if self.out_of_fuel() == 1 and (not self.has_separated):
            self.stage_separation()

        if self.out_of_fuel() == 2:
            self.stage.thrust = 0

        if self.altitutde_too_high():
            self.too_high = True
            return

        if track_state:
            self.state_trajectory.append(state)
            self.extra_trajectory.append(state_extra)
        
        return