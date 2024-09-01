"""
Define key global variables used throughout the other scripts
"""
import numpy as np

# Define global variables
G0 = 9.81 #m/s
MU = 3.986e14 #m^3/s^2
R_EARTH = 6350e3 #m

# Time values
DT = 1
T_FINAL = 10000

angle_inclination = np.deg2rad(10) #rad
angle_azimuth = np.deg2rad(0) #rad
t_start = 60 #s
t_burn = 10 #s

X0_PITCHOVER = [angle_inclination, angle_azimuth, t_start, t_burn]

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

def normalize(v):
    return v / np.linalg.norm(v)

def find_dcm(vec1, vec2):
    vec1 = normalize(vec1)
    vec2 = normalize(vec2)

    v = np.cross(vec1, vec2)
    s = np.linalg.norm(v)
    c = np.dot(vec1, vec2)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    if np.abs(s) < 1e-3:
        return np.eye(3)

    dcm = np.eye(3) + vx + vx.dot(vx) * ((1-c)/(s**2))

    return dcm