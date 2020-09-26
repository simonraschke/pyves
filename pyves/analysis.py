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



import pandas as pd
import numpy as np
from numpy.linalg import svd
from scipy import ndimage
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import normalize



def fullAnalysis(
        df : pd.DataFrame, 
        prms : dict,
        system = None
    ):
    assert("x" in df.columns)
    assert("y" in df.columns)
    assert("z" in df.columns)
    assert("ux" in df.columns)
    assert("uy" in df.columns)
    assert("uz" in df.columns)
    assert("sigma" in df.columns)
    assert("epsilon" in df.columns)
    assert("kappa" in df.columns)
    assert("gamma" in df.columns)

    analysis_prms = prms.get("analysis")
    DBSCAN_eps = analysis_prms.get("DBSCAN_eps", 1.2)
    neighbor_cutoff = analysis_prms.get("neighbor_cutoff", 1.2)
    volume_max_grid = analysis_prms.get("volume_max_grid", 80)
    volume_points_per_sigma = analysis_prms.get("volume_points_per_sigma", 3)

    system_prms = prms.get("system")
    box_prms = system_prms.get("box")
    dimensions = [box_prms["x"], box_prms["y"], box_prms["z"], 90, 90, 90]

    distances = distance_array(df.filter(["x","y","z"]).values, df.filter(["x","y","z"]).values, box=dimensions)
    sigma_scaling_mask = np.add.outer(df["sigma"].values, df["sigma"].values)/2

    
    df["cluster"] = DBSCAN(min_samples=2, eps=DBSCAN_eps, metric="precomputed", n_jobs=-1).fit(distances).labels_.astype(np.int16)
    
    unique, counts = np.unique(df["cluster"], return_counts=True)
    df["clustersize"] = df["cluster"].apply( lambda x: counts[np.where(unique == x)][0] if x >= 0 else 1 ).astype(np.uint16)

    df["subcluster"] = df.groupby("cluster").apply(lambda group: _subclusterLabels(group.name, group, DBSCAN_eps))["subcluster"].astype(np.uint8)

    df["neighbors"] = np.count_nonzero(distances <= sigma_scaling_mask*neighbor_cutoff, axis=1).astype(np.uint8) - 1

    df[["shiftx","shifty","shiftz"]] = df.groupby("cluster").apply(lambda group: shiftedCoordinates(group.name, group, dimensions))[["shiftx","shifty","shiftz"]].astype(np.float32)

    df["order"] = df.groupby("cluster").apply(lambda group: order(group.name, group))["order"].astype(np.float32)

    df["volume"] = df.groupby("cluster").apply(lambda group: volume(group.name, group, volume_max_grid, DBSCAN_eps, volume_points_per_sigma))["volume"].astype(np.float32)

    df["clustervolume"] = df[["clustersize","volume"]].prod(axis="columns")


    if isinstance(system, type(None)):
        system = makeSystem(df,prms)
        system.makeNeighborLists()
    df["epot"] = system.particleEnergies()
    df["chi"] = system.particleChiValues()

    return df



def makeSystem(
        df : pd.DataFrame, 
        prms : dict
    ):
    assert("x" in df.columns)
    assert("y" in df.columns)
    assert("z" in df.columns)
    assert("ux" in df.columns)
    assert("uy" in df.columns)
    assert("uz" in df.columns)
    import _pyves

    system = _pyves.System()
    system.box.x = prms["system"]["box"]["x"]
    system.box.y = prms["system"]["box"]["y"]
    system.box.z = prms["system"]["box"]["z"]
    system.interaction_cutoff = prms["system"]["interaction"]["cutoff"]
    system.cell_update_interval = prms["control"]["cell_update_interval"]
    system.neighbor_update_interval = prms["control"]["neighbor_update_interval"]
    system.neighbor_cutoff = prms["control"]["neighbor_cutoff"]

    for _, row in df.iterrows():
        system.particles.append(_pyves.Particle([row["x"], row["y"], row["z"]], [row["ux"], row["uy"], row["uz"]], 
            sigma=row["sigma"], kappa=row["kappa"], eps=row["epsilon"], gamma=row["gamma"], name=row["name"]))
        # system.particles[-1].initial_position = [row["initial_x"], row["initial_y"], row["initial_z"]]
        # system.particles[-1].initial_orientation = [row["initial_ux"], row["initial_uy"], row["initial_uz"]]
        # system.particles[-1].translation_bound_sq = row["translation_bound_sq"]
        # system.particles[-1].rotation_bound = row["rotation_bound"]
    
    # Setup cells
    box_dims = np.array([system.box.x, system.box.y, system.box.z])
    cells_per_dim = np.array(box_dims/prms["control"]["cell_min_size"]).astype(int)
    cell_actual_size = box_dims/cells_per_dim

    for x in np.arange(0, box_dims[0], cell_actual_size[0]):
        for y in np.arange(0, box_dims[1], cell_actual_size[1]):
            for z in np.arange(0, box_dims[2], cell_actual_size[2]):
                _min = np.array([x,y,z])
                _max = np.array([x,y,z]) + cell_actual_size
                system.cells.append(_pyves.Cell(_min, _max, box=system.box))

    for ci in system.cells:
        for cj in system.cells:
            if ci == cj:
                ci.region.append(cj)
            if ci.isNeighbourOf(cj):
                ci.proximity.append(cj)
                ci.region.append(cj)
    
    for c in system.cells:
        assert(c.assertIntegrity())

    cell_place_counter = 0
    for i, particle in enumerate(system.particles):
        for j, cell in enumerate(system.cells):
            if cell.insideCellBounds(particle):
                cell.particles.append(particle)
                # print(i, "to", j)
                cell_place_counter += 1
                break
    assert(cell_place_counter == len(system.particles))
    # assert(system.assertIntegrity())
    return system



def _subclusterLabels(label, group, eps):
    group["subcluster"] = 0
    if label == -1:
        return group
    else:
        # arange a DBSCAN without PBC to get subclusters
        coms = group.filter(["x","y","z"])
        distances_array = distance_array(coms.values, coms.values, box=None)
        group["subcluster"] = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=1).fit(distances_array).labels_
        return group



def shiftedCoordinates(label, group, dimensions):
    group[["shiftx","shifty","shiftz"]] = group[["x","y","z"]] 
    if label == -1:
        return group
    # get largest subcluster
    unique, counts = np.unique(group["subcluster"], return_counts=True)
    if len(unique) == 1:
        return group
    max_subclusterID = unique[counts == np.max(counts)][0]
    # calculate shifts per subcluster
    centers = group.groupby("subcluster")[['x','y','z']].mean()
    shifts = np.round(( -centers + centers.loc[max_subclusterID] )/dimensions[:3]).astype(int)
    shifts *= dimensions[:3]
    group[["shiftx","shifty","shiftz"]] += shifts.loc[group["subcluster"]].values
    # print(group)
    return group



def order(label, group):
    group["order"] = np.nan
    if label== -1:
        return group
    else:
        if group["gamma"].median() < 1e-3:
            if len(group.index) <= 3:
                return group
            shifted_coms = group.filter(['shiftx','shifty','shiftz'])
            normalized_orientations = group.filter(['ux','uy','uz'])
            center, normal = planeFit(shifted_coms.values.T)
            data = np.sum(normalized_orientations.values*normal, axis=1)
            if data.mean() < 0:
                data_compare = data * (-1)
                group["order"] = [data, data_compare][np.argmax([data.mean(), data_compare.mean()])]
            else:
                group["order"] = data
            return group
        else :
            shifted_coms = group.filter(['shiftx','shifty','shiftz'])
            normalized_orientations = group.filter(['ux','uy','uz'])
            data = np.sum(normalized_orientations.values*normalize(shifted_coms.sub(shifted_coms.mean()).values, copy=False), axis=1)
            group["order"] = data
            return group



def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]



def volume(label, group, max_points, eps, pps):
    """
    ppp : points per sigma
    eps : distance from point to count as in-volume
    """
    # volume = 0
    group["volume"] = np.pi * ((group["sigma"]*eps)**3) * 4/3
    if label == -1:
        return group
    else:
        xmin, xmax = group["shiftx"].min() - group["sigma"].max()*eps, group["shiftx"].max() + group["sigma"].max()*eps
        ymin, ymax = group["shifty"].min() - group["sigma"].max()*eps, group["shifty"].max() + group["sigma"].max()*eps
        zmin, zmax = group["shiftz"].min() - group["sigma"].max()*eps, group["shiftz"].max() + group["sigma"].max()*eps

        x_vector = np.linspace(xmin, xmax, np.min([int((xmax-xmin)*pps), max_points]), dtype=np.float16)
        y_vector = np.linspace(ymin, ymax, np.min([int((ymax-ymin)*pps), max_points]), dtype=np.float16)
        z_vector = np.linspace(zmin, zmax, np.min([int((zmax-zmin)*pps), max_points]), dtype=np.float16)

        gridpoints = np.vstack(np.meshgrid(x_vector, y_vector, z_vector)).reshape(3,-1).T
        #calculate the distance array with centres of masses of particles
        coms_cluster = group.filter(['shiftx','shifty','shiftz'])
        volume.distances_array_volume = distance_array(gridpoints, coms_cluster.values, box=None).astype(np.float32)
        mask = np.repeat([group["sigma"].values*eps], gridpoints.shape[0], axis=0)
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(volume.distances_array_volume <= mask, True, False).any(axis=1).reshape(x_vector.shape[0], y_vector.shape[0], z_vector.shape[0])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        group["volume"] = np.diff(x_vector)[0] * np.diff(y_vector)[0] * np.diff(z_vector)[0] * np.count_nonzero(isclose) / group["clustersize"]
        return group