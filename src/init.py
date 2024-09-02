from rocket import Rocket, RocketStage, PitchOver
import numpy as np
from globals import X0_PITCHOVER

EXAMPLE_PITCHOVER = PitchOver(X0_PITCHOVER)

def create_rocket(pitchover=EXAMPLE_PITCHOVER):
    # Initialize state
    # Payload
    mass_payload = 19000 #kg

    # Stage 1
    mass_struct = 25600 #kg
    mass_prop = 395700 #kg
    mass = mass_struct + mass_prop #kg
    Isp = 282 #s
    thrust = 7600e3 #N
    stage1 = RocketStage(mass=mass,
                         mass_struct=mass_struct,
                         mass_prop=mass_prop,
                         Isp=Isp,
                         thrust=thrust
                         )
    
    # Stage 2
    mass_struct = 3900 #kg
    mass_struct += mass_payload    
    mass_prop = 92670 #kg
    mass = mass_struct + mass_prop #kg
    Isp = 348 #s
    thrust = 981e3 #N
    stage2 = RocketStage(mass=mass,
                         mass_struct=mass_struct,
                         mass_prop=mass_prop,
                         Isp=Isp,
                         thrust=thrust
                         )
    
    # Pitchover is the design variable

    # Assume we launch out of Cape Canaveral (at midnight EST)
    v_launchpad = np.array([222, 70, -402])
    rx, ry, rz = 6377435, -1078282, 1329213 #m

    # v_launchpad = np.array([0, 0, 0])
    # rx, ry, rz = 0, 0, 6350e3
    vx, vy, vz = 0, 0, 0 #m/s
    mass = stage1.mass + stage2.mass #kg
    cd = 0.75
    area = np.pi/4 * 3.7**2

    init_state = np.array([rx, ry, rz, vx, vy, vz, mass])

    rocket = Rocket(init_state=init_state,
                    stage1=stage1,
                    stage2=stage2,
                    pitchover=pitchover,
                    cd=cd,
                    area=area,
                    v_launchpad=v_launchpad)
    
    return rocket