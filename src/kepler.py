# Filename: kepler.py
# Created by Koitchyy

import numpy as np
import coordinates as coord # use for coordinate transformations 

pi      = np.pi
deg2rad = pi / 180.0

def mean2TrueAnomaly(a, e, M0_rad, t0, t, mu):
    # Mean motion
    n = np.sqrt(mu / a ** 3)
    M  = (M0_rad + n * (t - t0)) % (2*pi) # wrap to 2π

    # Newton-Raphson method
    # Find eccentric anomaly from M = E - esin(E) 
    E = M.copy() # initial "guess"
    for _ in range(15):         
        tol = 1e-10
        # Newton step: ΔE = −f/f′
        deltaE = -(E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E     += deltaE
        # *** step-size convergence test ***
        if np.all(np.abs(deltaE) < tol):
            break

    # Use E (eccentric anomaly) to find nu (true anomaly)
    sin_nu =  np.sqrt(1-e**2) * np.sin(E) / (1 - e*np.cos(E))
    cos_nu = (np.cos(E) - e)      / (1 - e*np.cos(E))
    nu     = np.arctan2(sin_nu, cos_nu)
    return nu


def KOE2PQW(a, e, nu, mu):
    # Given the Keplerian Orbital Elements, will return the r-vector 
    # and velocity vector in the perifocal frame.
    r = a*(1 - e**2) / (1 + e * np.cos(nu)) # Found in problem 1a)
    x_PQW = r * np.cos(nu)
    y_PQW = r * np.sin(nu)
    z_PQW = 0.0
    v_x_PQW = -np.sqrt(mu/(a*(1-e**2))) * np.sin(nu)
    v_y_PQW = np.sqrt(mu/(a*(1-e**2))) * (e + np.cos(nu))
    v_z_PQW = 0.0

    # z components are zero, since r and v happen in orbital plane
    return (np.array([x_PQW, y_PQW, z_PQW]), 
            np.array([v_x_PQW, v_y_PQW, v_z_PQW]))

def elem2coord(a, e, i, omega, Omega, M0, t0, t, mu):
    #% INPUTS
    #% a - semimajor axis
    #% e - eccentricity
    #% i - inclination
    #% omega - argument of pericenter 
    #% Omega - longitude of the ascending node
    #% M0 - mean anomaly at time t0
    #% t0 - initial time
    #% t - time
    #% mu - gravitational constant (mu = GM)

    #% OUTPUT
    #% x,y,z - position
    #% vx, vy, vz – velocity
    ##########################

    i_rad     = i * deg2rad
    omega_rad = omega * deg2rad
    Omega_rad = Omega * deg2rad
    M0_rad    = M0 * deg2rad

    # Find the true anomaly
    nu = mean2TrueAnomaly(a, e, M0_rad, t0, t, mu)

    # Compute perifocal (PQW) r, v coords
    r_PQW, v_PQW = KOE2PQW(a, e, nu, mu)

    # Use coordinate transformations to convert to r and v
    R = (coord.rotation(-Omega_rad, 3) @
         coord.rotation(-i_rad,     1) @
         coord.rotation(-omega_rad, 3))

    # transpose to go PQW ➜ ECI
    x, y, z = R @ r_PQW
    vx, vy, vz = R @ v_PQW

    ##########################
    return x, y, z, vx, vy, vz

def coord2elem(x, y, z, vx, vy, vz, mu):

    r = np.array([x, y, z], dtype=float)
    v = np.array([vx, vy, vz], dtype=float)
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    h = np.array([y * vz - z * vy, 
                  z * vx - x * vz,
                  x * vy - y * vx])
    
    h_norm = np.linalg.norm(h)

    W_x = h[0] / h_norm
    W_y = h[1] / h_norm
    W_z = h[2] / h_norm

    # 2.58 from Gill and Montenbruck
    i = np.arctan2(np.sqrt(W_x**2 + W_y**2), W_z)
    Omega = np.arctan2(W_x, -W_y)

    # semi latus rectum
    p = h_norm ** 2 / mu

    # vis-viva law
    a = 1.0 / ((2.0 /r_norm) - (v_norm ** 2 / mu))

    # Mean motion
    n = np.sqrt(mu / a**3)

    # Eccentricity (a is always positive)
    e = np.sqrt(1 - p/a)

    # Find eccentric anomaly
    E = np.arctan2((r @ v)/((a*a)*n), (1 - r_norm/a)) % (2*pi)

    # Find Mean anomaly
    M = (E - e * np.sin(E)) % (2*pi)

    # Find argument of latitude u = nu + Ω
    u = np.arctan2(z, -x * W_y + y * W_x)

    # Find true anomaly nu
    nu = np.arctan2(np.sqrt(1- e**2) * np.sin(E), np.cos(E) - e)

    # Get argument of perigee
    omega = (u - nu) % (2*pi)

    return (a, e, i, omega, Omega, M)