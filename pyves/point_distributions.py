# pyves - python3 bindings for an easy use of vesicle2
# Copyright (C) 2020 Simon Raschke

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



import numpy as np


def fibonacci_sphere_points(samples=1):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append([x, y, z])

    return np.array(points)



def sunflower_sphere_points(samples=1):
    indices = np.arange(0, samples, dtype=float) + 0.5

    phi = np.arccos(1 - 2*indices/samples)
    theta = np.pi * (1 + 5**0.5) * indices
    
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi);
    return np.array([x,y,z]).transpose()



def grid_plane_points(samples=1):
    per_dim = int(np.sqrt(samples))
    if per_dim**2 != samples:
        raise RuntimeError(f"Got {samples} samples, but need int**2 samples.")
    xvals = np.linspace(-1, 1, per_dim, endpoint=True)
    yvals = np.linspace(-1, 1, per_dim, endpoint=True)
    zvals = np.array([0])
    grid = np.meshgrid(xvals, yvals, zvals)
    points = np.array(grid).transpose().reshape(per_dim**2, 3)
    return points



def hexagonal_lattice_points(x=1.0, y=1.0, distance=1.0):
    def rotate(vec, angle_in_deg):
        from scipy.spatial.transform import Rotation as R
        return R.from_euler('z', angle_in_deg, degrees=True).apply(vec)
        
    baserow = np.array([np.arange(0, x, distance), np.zeros_like(np.arange(0, x, distance)), np.zeros_like(np.arange(0, x, distance))]).T
    rowheight = rotate([distance,0,0], 30)[0]
    rowoffset = rotate([distance,0,0], 60)[0]
    y_heights = np.arange(0, y, rowheight)
    grid = baserow.copy()
    for i, y_val in enumerate(y_heights[1:]):
        newrow = baserow.copy()
        newrow[:,1] = y_val
        if i % 2 == 0:
            newrow[:,0] += rowoffset
        grid = np.append(grid, newrow, axis=0)    
    return grid