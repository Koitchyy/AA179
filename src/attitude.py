import numpy as np
import coordinates as coord

def inertia_cube(m, a):
    I = 1/6*m*a**2*np.array([[1, 0, 0], 
                             [0, 1, 0],
                             [0, 0, 1]])
    return I

def inertia_plate(m, a, b):
    I = 1/12*m*np.array([[b**2, 0,    0], 
                         [0,    a**2, 0],
                         [0,    0,    a**2+b**2]])
    return I

def inertia_box(m, a, b, c):
    I = 1/12*m*np.array([[b**2+c**2, 0,         0], 
                         [0,         a**2+c**2, 0],
                         [0,         0,         a**2+b**2]])
    return I

def inertia_translate(m, I, dx, dy, dz):
    # this is for the parallel axis theorem
    dI = m*np.array([[dy**2+dz**2, -dx*dy,      -dx*dz], 
                     [-dx*dy,      dx**2+dz**2, -dy*dz],
                     [-dx*dz,      -dy*dz,      dx**2+dy**2]])
    return I+dI

def principal_inertia(I):
    # diagonalize moment of inertia tensor here
    I_prin, axes = np.linalg.eig(I)
    return I_prin, axes

def inertia_rotate(I, theta, ax):
    # compute moment of inertia tensor in a rotated frame
    # rotate around axis (1 for x, 2 for y, 3 for z) by angle theta
    R = coord.rotation(theta, ax)
    I_rot = R @ I @ np.transpose(R)
    return I_rot
