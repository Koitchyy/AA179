#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 16:50:18 2024

@author: eai
"""

import numpy as np

def elem2coord(a,e,i,w,W,M0,t0,t,mu):

    #% INPUTS
    #% a - semimajor axis
    #% e - eccentricity
    #% i - inclination (radians)
    #% w - argument of pericenter (pericenter is the same as periapsis)
    #% W - longitude of the accending nodes
    #% M0 - mean anomaly at time t0
    #% t0 - initial time
    #% t  - time
    #% mu  - gravitational constant (GM)

    #% OUTPUT
    #% x,y,z      - position
    #% vx, vy, vz - velocity
    
    n = np.sqrt(mu / a / a / a)
    M = M0 + n * (t - t0)

    # solving Kepler equation
    d = np.array(1.0)
    E = M + e * np.sin(M)
    eps = 1e-10
    while (d > eps).any():
        E1 = E
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        d = abs(E - E1)

    v = np.sqrt((1 + e) / (1 - e)) 
    v = 2 * np.arctan(v * np.tan(E / 2))

    u = w + v

    ax = np.cos(W) * np.cos(u) - np.sin(W) * np.sin(u) * np.cos(i)
    ay = np.sin(W) * np.cos(u) + np.cos(W) * np.sin(u) * np.cos(i)
    az = np.sin(u) * np.sin(i)
    
    r = a * (1 - e * np.cos(E))

    x = r * ax 
    y = r * ay 
    z = r * az

    axu = -np.cos(W) * np.sin(u) - np.sin(W) * np.cos(u) * np.cos(i)
    ayu = -np.sin(W) * np.sin(u) + np.cos(W) * np.cos(u) * np.cos(i)
    azu =  np.cos(u) * np.sin(i)

    p = a * (1 - e * e)

    # velocities
    vx = np.sqrt(mu / p) * ((e * np.sin(v)) * ax + (1 + e * np.cos(v)) * axu)
    vy = np.sqrt(mu / p) * ((e * np.sin(v)) * ay + (1 + e * np.cos(v)) * ayu)
    vz = np.sqrt(mu / p) * ((e * np.sin(v)) * az + (1 + e * np.cos(v)) * azu)
    
    return x, y, z, vx, vy, vz

def coord2elem(x, y, z, vx, vy, vz, mu):
    
    r = np.sqrt(x  **2  + y  ** 2 + z ** 2)
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    
    h = (v ** 2) / 2 - mu / r
    
    c1 = y * vz - z * vy
    c2 = z * vx - x * vz
    c3 = x * vy - y * vx
    
    l1 = - mu * x / r + vy * c3 - vz * c2
    l2 = - mu * y / r + vz * c1 - vx * c3
    l3 = - mu * z / r + vx * c2 - vy * c1
    
    c = np.sqrt(c1 * c1 + c2 * c2 + c3 * c3)
    l = np.sqrt(l1 * l1 + l2 * l2 + l3 * l3)
    
    a = - mu / 2.0 / h
    e = l / mu
    
    # getting inclination
    cosi = c3 / c
    sini = np.sqrt(1 - cosi ** 2)
    i    = np.arctan2(sini, cosi)
    
    # getting longitude of the ascending node
    W = np.arctan2(c1 / c /sini, -c2 / c / sini)
    
    # getting argument of the pericenter
    w = np.arctan2(l3 / l / sini, l1 / l * np.cos(W) + l2 / l * np.sin(W))
    
    # getting u
    u = np.arctan2(z / r / sini, x / r * np.cos(W) + y / r * np.sin(W))
    
    # getting true anomaly
    v = np.arctan2(np.sin(u) * np.cos(w) - np.cos(u) * np.sin(w), np.cos(u) * np.cos(w) + np.sin(u) * np.sin(w))
    
    # getting eccentric anomaly
    sinE = np.sqrt(1-e*e) * np.sin(v) / (1 + e * np.cos(v))
    cosE = (np.cos(v) + e) / (1 + e * np.cos(v))
    E = v + np.arctan((sinE * np.cos(v) - cosE * np.sin(v)) / (cosE * np.cos(v) + sinE * np.sin(v)))
    
    # getting mean anomaly
    M = E - e * np.sin(E)
    
    return (a, e, i, np.mod(w, 2*np.pi), np.mod(W, 2*np.pi), np.mod(M, 2*np.pi))

    



    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    