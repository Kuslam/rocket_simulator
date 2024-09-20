import numpy as np
from globals import *

def get_orbital_elements(r, v):
    """
    Determine the following orbital elements based on a given position (r) and velocity vector (v):
    - a: apoapsis
    - e: eccentricity
    - i: inclination
    - Om: Right Ascension of the Ascending Node (RAAN)
    - w: argument of periapsis
    - nu: true anamoly
    """

    # Find norms
    v_norm = np.linalg.norm(v)
    r_norm = np.linalg.norm(r)

    # Find key vectors
    h = np.cross(r,v)
    n = np.cross(np.array([0,0,1]),h)
    e_vec = ((v_norm**2 - MU/r_norm)*r - np.dot(r,v)*v)/MU
    e = np.linalg.norm(e_vec)
    
    # Find energy in orbit
    E = v_norm**2/2 - MU/r_norm

    # Find a and p
    if e <= 1:
        a = -MU/(2*E)
        p = a * (1-e**2)
    else:
        a = np.inf
        p = np.linalg.norm(h)**2/MU

    i = np.arccos(h[2]/np.linalg.norm(h))
    Om = np.arccos(n[0]/np.linalg.norm(n))
    w = np.arccos(np.dot(n,e_vec)/(np.linalg.norm(n)*e))
    nu = np.arccos(np.dot(r,e_vec)/(np.linalg.norm(r)*e))

    # Final checks
    if n[1] < 0:
        Om = 2*np.pi - Om
    if e_vec[1]:
        w = 2*np.pi - w
    if np.dot(r,v) < 0:
        nu = 2*np.pi - nu

    return a, e, i, Om, w, nu

