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
    # grid = np.meshgrid(xx, yy, zz)
    points = np.array(grid).transpose().reshape(per_dim**2, 3)
    return points
