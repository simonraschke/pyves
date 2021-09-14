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



import os
import re
import sys
import time
import json
import tables
import pandas as pd
import numpy as np

from pyves import utility
from .signal_handler import SignalHandler, ProgramState

from concurrent.futures import ProcessPoolExecutor, as_completed
from numpy.linalg import svd
from scipy import ndimage
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import normalize



def analyzeTrajectory(
    inpath = None,
    outpath = None,
    prmspath = "parameters.json",
    threads = -1,
    timestats = False
):
    if not isinstance(inpath, type(None)):
        inpath = os.path.abspath(inpath)
    else:
        inpath = ""
    if not isinstance(outpath, type(None)):
        outpath = os.path.abspath(outpath)
    else:
        outpath = ""
        
    prmspath = os.path.abspath(prmspath)
    with open(prmspath, 'r') as prms_file:
        prms = json.loads(prms_file.read())
    if threads == -1:
        threads = prms["hardware"].get("threads", 1)

    if inpath == outpath and len(outpath) > 0:
        pass
    elif len(outpath) == 0 and os.path.exists(inpath):
        outpath = inpath
    elif len(inpath) == 0 and len(outpath) > 0:
        inpath = os.path.join(prms["output"]["dir"], prms["output"]["filename"])
    elif len(outpath) == 0 and len(inpath) == 0:
        inpath = os.path.join(prms["output"]["dir"], prms["output"]["filename"])
        outpath = os.path.join(prms["output"]["dir"], prms["output"]["filename"])
    elif os.path.exists(inpath) and not os.path.exists(outpath):
        pass
    else:
        raise NotImplementedError(f"{inpath}, {outpath}")

    print("analyze trajectory")
    print("input      path:", inpath)
    print("output     path:", outpath)
    print("parameters path:", prmspath)
    
    tables.file._open_files.close_all()
    store = pd.HDFStore(inpath, mode="r")
    keys = [s for s in store.keys() if s.startswith("/time")]
    keys = sorted(keys, key=lambda x:int(re.findall('(?<=time)\d+', x)[0]))
    store.close()
    store = None


    for group in utility.chunker(keys, threads):
        # dfchunk is (key, df, metadata)
        dfchunk = [(key, *utility.h5load(inpath, key=key)) for key in group]
        with ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_key = {}
            for key, df, metadata in dfchunk:
                assert isinstance(metadata, dict)
                if not SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    future_to_key[executor.submit(analyzeSnapshot, df=df, prms=prms, metadata=metadata, timestats=True)] = (key, metadata)
                else:
                    break
            for future in as_completed(future_to_key):
                if not SignalHandler.ProgramState is ProgramState.SHUTDOWN:
                    key, metadata = future_to_key[future]
                    assert isinstance(metadata, dict)
                    # try:
                    #     df, execution_time = future.result()
                    #     if timestats: print(f"analysis took {execution_time:.2f} s | written to {outpath}{key}")
                    #     utility.h5store(outpath, key, df, **metadata)
                    # except Exception as e:
                    #     print(e)
                    df, execution_time = future.result()
                    if timestats: print(f"analysis took {execution_time:.2f} s | written to {outpath}{key}")
                    utility.h5store(outpath, key, df, **metadata)
                else:
                    break



def analyzeSnapshot(
        df : pd.DataFrame, 
        prms : dict,
        metadata : dict,
        system = None,
        timestats = False
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
    assert isinstance(prms, dict)
    assert isinstance(metadata, dict)
    
    starttime = time.perf_counter()

    analysis_prms = prms.get("analysis")    
    DBSCAN_eps = analysis_prms.get("DBSCAN_eps", 1.2)
    neighbor_cutoff = analysis_prms.get("neighbor_cutoff", 1.2)
    curvature_cutoff = analysis_prms.get("curvature_cutoff", 1.2)
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

    df["curvature"] = curvature(df, dimensions, cutoff=curvature_cutoff)


    if isinstance(system, type(None)):
        system = utility._makeSystem(df, prms, metadata)
        system.makeNeighborLists()
    df["epot"] = system.particleEnergies()
    df["chi"] = system.particleChiValues()
    df["surfacepot"] = system.particleSurfacePotentialValues()
    df["externalpot"] = system.particleExternalPotentialValues()
    try:
        assert df[(df.z < 3) & (df.surface_affinity_translation > 1e-3)]['surfacepot'].ge(1e5).all()
    except AssertionError as e:
        print(df[(df.z < 3) & (df.surface_affinity_translation > 1e-3)][["z","surfacepot"]])
        raise e

    endtime = time.perf_counter()
    if timestats: 
        return df, endtime-starttime
    else:
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



def _volume_calculation(shiftx, shifty, shiftz, sigma, clustersize, max_points, eps, pps):
    xmin = np.min(shiftx) - np.max(shiftx)*eps 
    xmax = np.max(shiftx) + np.max(shiftx)*eps
    ymin = np.min(shifty) - np.max(shifty)*eps 
    ymax = np.max(shifty) + np.max(shifty)*eps
    zmin = np.min(shiftz) - np.max(shiftz)*eps 
    zmax = np.max(shiftz) + np.max(shiftz)*eps

    x_vector = np.linspace(xmin, xmax, np.min(np.array([int((xmax-xmin)*pps), max_points])))
    y_vector = np.linspace(ymin, ymax, np.min(np.array([int((ymax-ymin)*pps), max_points])))
    z_vector = np.linspace(zmin, zmax, np.min(np.array([int((zmax-zmin)*pps), max_points])))

    gridpoints = np.vstack(np.meshgrid(x_vector, y_vector, z_vector)).reshape(3,-1).T
    #calculate the distance array with centres of masses of particles
    # coms_cluster = group.filter(['shiftx','shifty','shiftz'])
    coms_cluster = np.row_stack([shiftx,shifty,shiftz]).T
    # print(coms_cluster)
    distances_array_volume = distance_array(gridpoints, coms_cluster, box=None)

    mask = np.repeat([sigma*eps], gridpoints.shape[0], axis=0)
    # check if any point in distance array row is close enough, then reshape to meshgrid
    # result is a binary meshgrid with 1 for the cluster shell region
    isclose = np.where(distances_array_volume <= mask, True, False).any(axis=1).reshape(x_vector.shape[0], y_vector.shape[0], z_vector.shape[0])
    # fill hole inside the shell region
    isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
    # calc volum from all points inside cluster
    return np.diff(x_vector)[0] * np.diff(y_vector)[0] * np.diff(z_vector)[0] * np.count_nonzero(isclose) / clustersize



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
        distances_array_volume = distance_array(gridpoints, coms_cluster.values, box=None).astype(np.float32)
        mask = np.repeat([group["sigma"].values*eps], gridpoints.shape[0], axis=0)
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(distances_array_volume <= mask, True, False).any(axis=1).reshape(x_vector.shape[0], y_vector.shape[0], z_vector.shape[0])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        group["volume"] = np.diff(x_vector)[0] * np.diff(y_vector)[0] * np.diff(z_vector)[0] * np.count_nonzero(isclose) / group["clustersize"]
        return group
        # coms_cluster = np.array([group['shiftx'].values, group['shifty'].values, group['shiftz'].values])

        # group["volume"] = _volume_calculation(
        #     group['shiftx'].values, 
        #     group['shifty'].values, 
        #     group['shiftz'].values, 
        #     group['sigma'].values,
        #     group['clustersize'].values, 
        #     max_points, eps, pps
        # )
        # return group



def getPairs(distances_array, max_cutoff, same=False):
    valid = distances_array < max_cutoff
    np.fill_diagonal(valid, same)
    pairs = np.column_stack(np.where(valid))
    if len(pairs) == 0:
        return []
    pairs = np.sort(pairs, axis=1)
    _, index = np.unique(pairs, axis=0, return_index=True)
    return pairs[index]



def curvature(particledata, dimensions, cutoff):
    coms = particledata.filter(['shiftx','shifty','shiftz'])
    orientations = particledata.filter(['ux','uy','uz'])
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    pairs = getPairs(distances_array, cutoff)
    if len(pairs) == 0:
        return np.nan
    origin_orientations = orientations.values[pairs[:,0]]
    origin_connections = np.subtract(coms.values[pairs[:,1]], coms.values[pairs[:,0]])
    # origin_connections = np.divide(origin_connections, 10)
    projections = np.einsum('ij,ij->i', origin_connections, origin_orientations) # same as (origin_connections * origin_orientations).sum(axis=1) BUT FASTER
    projections_array = np.zeros_like(distances_array)
    pairs_t = pairs.T
    projections_array[tuple(pairs_t)] = projections
    projections_array[tuple([pairs_t[1], pairs_t[0]])] = projections
    sums = np.sum(projections_array, axis=1)
    nums = np.count_nonzero(projections_array, axis=1)
    averages = np.zeros_like(sums)
    averages[np.where(nums>0)] = sums[np.where(nums>0)]/nums[np.where(nums>0)]
    
    coms = particledata.filter(['x','y','z'])
    distances_array = distance_array(coms.values, coms.values, box=None)
    pairs = getPairs(distances_array, cutoff)
    if len(pairs) == 0:
        return np.full((len(particledata.index),1), np.nan, dtype=np.float32)
    origin_orientations = orientations.values[pairs[:,0]]
    origin_connections = np.subtract(coms.values[pairs[:,1]], coms.values[pairs[:,0]])
    # origin_connections = np.divide(origin_connections, 10)
    projections = np.einsum('ij,ij->i', origin_connections, origin_orientations) # same as (origin_connections * origin_orientations).sum(axis=1) BUT FASTER
    projections_array = np.zeros_like(distances_array)
    pairs_t = pairs.T
    projections_array[tuple(pairs_t)] = projections
    projections_array[tuple([pairs_t[1], pairs_t[0]])] = projections
    _sums = np.sum(projections_array, axis=1)
    _nums = np.count_nonzero(projections_array, axis=1)
    _averages = np.zeros_like(_sums)
    _averages[np.where(_nums>0)] = _sums[np.where(_nums>0)]/nums[np.where(_nums>0)]

    condition = np.where(np.logical_or(averages<-1, averages>1))
    
    averages[condition] = _averages[condition]
    return np.nan_to_num(averages)