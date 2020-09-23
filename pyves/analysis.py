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
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN



def fullAnalysis(
        df : pd.DataFrame, 
        prms : dict
    ):
    assert("x" in df.columns)
    assert("y" in df.columns)
    assert("z" in df.columns)

    analysis_prms = prms.get("analysis")
    DBSCAN_eps = analysis_prms["DBSCAN_eps"]
    neighbor_cutoff = analysis_prms["neighbor_cutoff"]

    system_prms = prms.get("system")
    box_prms = system_prms.get("box")
    dimensions = [box_prms["x"], box_prms["y"], box_prms["z"], 90, 90, 90]

    distances = distance_array(df.filter(["x","y","z"]).values, df.filter(["x","y","z"]).values, box=dimensions)
    
    df["cluster"] = DBSCAN(min_samples=2, eps=DBSCAN_eps, metric="precomputed", n_jobs=-1).fit(distances).labels_.astype(np.int32)
    df["subcluster"] = df.groupby("cluster").apply(lambda group: _subclusterLabels(group.name, group, DBSCAN_eps))["subcluster"].astype(np.uint8)
    df["neighbors"] = np.count_nonzero(distances <= neighbor_cutoff, axis=1).astype(np.int16) - 1
    # print("\n",df)

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
