#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:49:08 2024

@author: eai
"""
import numpy as np

def cart2sph(x,y,z):
    azimuth   = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    return azimuth, elevation, r

def sph2cart(azimuth,elevation,r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z


def BLH2XYZ(B, L, H, a, e):
    
    N = a / np.sqrt(1 - (e * np.sin(B)) ** 2)
    
    x = (N + H)             * np.cos(B) * np.cos(L)
    y = (N + H)             * np.cos(B) * np.sin(L)
    z = (N + H - N * e * e) * np.sin(B)
    return x, y, z

def XYZ2BLH(x, y, z, a, e):
    e2     = e * e
    r_perp = np.sqrt(x * x + y * y)
    L      = np.arctan2(y,x)
    B      = np.arctan(z / r_perp)
    tol = 1e-10
    eps = np.array(1.0)
    while (eps > tol).any():
        B_old = B
        B     = np.arctan((z / r_perp + (a * e2 * np.sin(B_old)) / r_perp / np.sqrt(1 - e2 * np.sin(B_old) ** 2)))
        eps = np.abs(B - B_old)

    H = r_perp * np.cos(B) + z * np.sin(B) - a * np.sqrt((1 - e2*np.sin(B)**2))
    
    return B, L, H

def rotation(theta, ax):
    if ax == 1:
        R = np.array([[1, 0,           0], \
                      [0, np.cos(theta),  np.sin(theta)],
                      [0, -np.sin(theta), np.cos(theta)]])
    elif ax == 2:
        R = np.array([[np.cos(theta), 0,           -np.sin(theta)], \
                      [0,             1,            0],
                      [np.sin(theta), 0, np.cos(theta)]])
   
    elif ax == 3:
        R = np.array([[np.cos(theta),  np.sin(theta),  0], \
                      [-np.sin(theta), np.cos(theta),  0],
                      [0,              0,              1]])
    else:
        print("wrong axis number")
        
    return R


def XYZ2NEU(x, y, z, xc, yc, zc, a, e):
    
    (B, L, H) = XYZ2BLH(xc, yc, zc, a, e)
    
    x_rel = x - xc
    y_rel = y - yc
    z_rel = z - zc
    
    # calculate NEU
    N = -np.sin(B) * np.cos(L) * x_rel -\
         np.sin(B) * np.sin(L) * y_rel +\
         np.cos(B) * z_rel
    
    E = -np.sin(L) * x_rel + \
         np.cos(L) * y_rel 
         
    U =  np.cos(B) * np.cos(L) * x_rel +\
         np.cos(B) * np.sin(L) * y_rel +\
         np.sin(B) * z_rel
         
    return N, E, U
    

def gregorian2JD(month, day, year):
    """Convert the given Gregorian date to the equivalent Julian Day Number."""
    if month < 1 or month > 12:
        raise ValueError("month must be between 1 and 12, inclusive")
    if day < 1 or day > 31:
        raise ValueError("day must be between 1 and 31, inclusive")
    A = int((month - 14) / 12)
    B = 1461 * (year + 4800 + A)
    C = 367 * (month - 2 - 12 * A)
    E = int((year + 4900 + A) / 100)
    JDN = int(B / 4) + int(C / 12) - int(3 * E / 4) + day - 32075.5
    return JDN
    

    
    
        
        
        
        