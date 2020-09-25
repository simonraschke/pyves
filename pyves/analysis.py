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
        prms : dict
    ):
    assert("x" in df.columns)
    assert("y" in df.columns)
    assert("z" in df.columns)

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
    df["subcluster"] = df.groupby("cluster").apply(lambda group: _subclusterLabels(group.name, group, DBSCAN_eps))["subcluster"].astype(np.uint8)
    df["neighbors"] = np.count_nonzero(distances <= sigma_scaling_mask*neighbor_cutoff, axis=1).astype(np.int16) - 1
    df[["shiftx","shifty","shiftz"]] = df.groupby("cluster").apply(lambda group: shiftedCoordinates(group.name, group, dimensions))[["shiftx","shifty","shiftz"]].astype(np.float32)
    df["order"] = df.groupby("cluster").apply(lambda group: order(group.name, group))["order"].astype(np.float32)
    # df["volume"] = df.groupby("cluster").apply(lambda group: volume(group.name, group, volume_max_grid, DBSCAN_eps, volume_points_per_sigma))["volume"].astype(np.float32)
    print("\n",df)
    print(df.info())

    return df



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

        xx, yy, zz = np.meshgrid(x_vector, y_vector, z_vector)
        meshgrid = np.meshgrid(x_vector, y_vector, z_vector)
        # stack them together as array of 3D points
        meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
        print(meshgrid.shape)
        #calculate the distance array with centres of masses of particles
        coms_cluster = group.filter(['shiftx','shifty','shiftz'])
        volume.distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None).astype(np.float32)
        # make sclaing for individual sigmas
        sigma_scaling_mask = np.add.outer(group["sigma"].values, group["sigma"].values)/2
        print("\n\ncluster size", group.index.size)
        print(np.unique(meshgrid[:,0]))
        print(np.unique(meshgrid[:,0]).shape)
        print(np.unique(meshgrid[:,1]).shape)
        print(np.unique(meshgrid[:,2]).shape)
        print("max_points", max_points)
        print("pps", pps)
        print("meshgrid", meshgrid.shape)
        print("coms_cluster", coms_cluster.shape)
        print("distances_array_volume", volume.distances_array_volume.shape)
        print("sigma_scaling_mask", sigma_scaling_mask.shape)
        print("\n")
        sys.exit()
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(volume.distances_array_volume <= sigma_scaling_mask*eps, True, False).any(axis=1).reshape(xx.shape[0], yy.shape[1], zz.shape[2])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        group["volume"] = ((1.0/pps)**3)*np.count_nonzero(isclose)
        return group