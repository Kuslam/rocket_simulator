from rocket import Rocket, RocketStage
import numpy as np

def create_rocket():
    # Initialize state
    # Payload
    mass_payload = 100000 #kg

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
    mass_prop = 92670 #kg
    mass = mass_struct + mass_prop + mass_payload #kg
    Isp = 348 #s
    thrust = 981e3 #N
    stage2 = RocketStage(mass=mass,
                         mass_struct=mass_struct,
                         mass_prop=mass_prop,
                         Isp=Isp,
                         thrust=thrust
                         )

    # Assume we launch out of Cape Canaveral
    v_launchpad = np.array([0, 0, 0])

    # rx, ry, rz = 916.361e3, -5539.907e3, 3014.812e3 #m
    rx, ry, rz = 0, 0, 6350e3
    vx, vy, vz = 0, 0, 0 #m/s
    mass = stage1.mass + stage2.mass #kg
    cd = 0.75
    area = np.pi/4 * 3.7**2

    init_state = np.array([rx, ry, rz, vx, vy, vz, mass])

    rocket = Rocket(init_state=init_state,
                    stage1=stage1,
                    stage2=stage2,
                    cd=cd,
                    area=area,
                    v_launchpad=v_launchpad)
    
    return rocket