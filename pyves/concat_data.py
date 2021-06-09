import os
import re
import gc
import time
import hashlib
import sqlite3
import concurrent.futures

import numpy as np
import pandas as pd
import datetime as dt

from pathlib import Path

from .utility import h5load, _h5load_inner



def _yield_states_from_HDF(
    path,
    threads = 1,
    key_prefix = "/time",
):
    store = pd.HDFStore(path, mode="r")
    keys = [s for s in store.keys() if s.startswith(key_prefix)]
    keys = sorted(keys, key=lambda x:int(re.findall(f"(?<={key_prefix})\d+", x)[0]))
    store.close()
    store = None

    for key in keys:
        yield h5load(path, key)



def _load_states_from_HDF(
    path,
    threads = 1,
    key_prefix = "/time",
):
    with pd.HDFStore(path, mode="r") as store:
        keys = [s for s in store.keys() if s.startswith(key_prefix)]
        keys = sorted(keys, key=lambda x:int(re.findall(f"(?<={key_prefix})\d+", x)[0]))

        data_meta = [(_h5load_inner(store, key)) for key in keys]
        
        # print(data_meta)
        return data_meta



def analyzeStates(
    path,
    threads = 1,
    key_prefix = "/time",
    sys = True,
    clstr = True,
    clstr_min_size = 30
):
    print("[ONGOING] ", path)

    store = pd.HDFStore(path, mode="r")
    keys = [s for s in store.keys() if s.startswith(key_prefix)]
    keys = sorted(keys, key=lambda x:int(re.findall(f"(?<={key_prefix})\d+", x)[0]))
    # print(keys)
    store.close()
    store = None

    if sys:
        sysdata = []
    if clstr:
        clstrdata = []
    
    path_str_for_hash = re.sub("[:]", "", str(path)).lstrip("CDEFGH").replace('\\', '/')

    p = hashlib.md5()
    p.update(path_str_for_hash.encode('utf-8'))
    path_id = p.hexdigest()

    for df, meta in _load_states_from_HDF(path, 1, key_prefix):

        df_free  = df[df["clustersize"].le(clstr_min_size)]
        df_clstr = df[~df["clustersize"].le(clstr_min_size)]
        assert(df_clstr.index.size + df_free.index.size == df.index.size)

        if sys:
            # print(df[["translation_bound_sq", "rotation_bound"]])
            m = hashlib.md5()
            m.update(path_str_for_hash.encode('utf-8'))
            m.update(str(meta["time"]).encode('utf-8'))
            sys_identifier = m.hexdigest()
            
            data = {}
            data["id"] = sys_identifier
            data["simulation"] = path_id
            
            data["time"] = meta["time"]
            data["gamma"] = df["gamma"].mean()
            data["sigma"] = df["sigma"].mean() 
            data["kappa"] = df["kappa"].mean()
            data["epsilon"] = df["epsilon"].mean()
            data["x"] = meta["box.x"]
            data["y"] = meta["box.y"]
            data["z"] = meta["box.z"]
            data["volume"] = data["x"]*data["y"]*data["z"]
            data["temperature"] = np.round(meta["temperature"], 3)
            data["cluster_min_size"] = clstr_min_size
            data["walltime"] = meta['simulation_walltime']
            data['translation_step'] = meta['translation_step']
            data['rotation_step'] = meta['rotation_step']
            data['interaction_cutoff' ] = meta['interaction.cutoff']
            data['interaction_surface'] = meta['interaction.surface']
            data['interaction_surface_width'] = meta['interaction.surface_width']

            data["size"] = df.index.size
            data["size_free"] = df_free.index.size
            data["size_restricted"] = df[(df["translation_bound_sq"] < 1e10) | (df["rotation_bound"] < 1e10)].index.size
            data["size_trans_restricted"] = df[df["translation_bound_sq"] < 1e10].index.size
            data["size_rot_restricted"] = df[df["rotation_bound"] < 1e10].index.size
            data["size_mobile"] = df[(df["translation_bound_sq"] >= 1e10) & (df["rotation_bound"] >= 1e10)].index.size
            data["size_cluster"] = df_clstr.index.size
            data["clustervolume_cumulated"] = df_clstr.groupby("cluster")["clustervolume"].mean().sum() if data["size_cluster"] > 0 else 0.0
            data["density"] = data["size"]/data["volume"]
            data["density_free"] = float(data["size_free"])/(data["volume"]-data["clustervolume_cumulated"])
            
            data["epot"] = df["epot"].mean()
            data["epot_cluster"] = df_clstr["epot"].mean()
            data["surfacepot"] = df["surfacepot"].mean()
            data["surfacepot_cluster"] = df_clstr["surfacepot"].mean()
            data["externalpot"] = df["externalpot"].mean()
            data["externalpot_cluster"] = df_clstr["externalpot"].mean()
            
            data["chi"] = df["chi"].mean()
            data["chi_cluster"] = df_clstr["chi"].mean()
            data["order_simple"] = df["order"].mean()
            data["order_cluster"] = df_clstr["order"].mean()
            data["neighbors"] = df["neighbors"].mean()
            data["neighbors_cluster"] = df_clstr["neighbors"].mean()
            data["clustersize"] = df_clstr["clustersize"].mean()
            data["clustervolume"] = df_clstr["clustervolume"].mean()

            if data['interaction_surface'] and data['interaction_surface_width']:
                surface_height = (meta["box.z"]-meta['interaction.surface_width'])
                df_surface = df[(df["z"] > surface_height)]
                df_nonsurface = df[(df["z"] <= surface_height)]
                df_nonsurface_free = df_nonsurface[df_nonsurface["clustersize"].le(clstr_min_size)]
                df_nonsurface_clstr = df_nonsurface[~df_nonsurface["clustersize"].le(clstr_min_size)]
                data["N_on_z_surface"] = df_surface.index.size
                data["area_pP_z_surface"] = df_surface["sigma"].mean() * (2.**(1./6)/2)**2 * (2.*np.sqrt(3))
                data["area_occupied_z_surface"] = data["area_pP_z_surface"] * data["N_on_z_surface"]
                data["max_coverage_z_surface"] = 1./data["area_pP_z_surface"]
                data["coverage_z_surface"] = data["area_occupied_z_surface"] / (data["x"]*data["y"]) / data["max_coverage_z_surface"]
                data["bulk_volume"] = data["x"] * data["y"] * (meta["cell_min_size"] + data['interaction_surface_width'])
                data["bulk_density"] = float(df_nonsurface.index.size) / data["bulk_volume"]
                data["bulk_clustervolume_cumulated"] = df_nonsurface_clstr.groupby("cluster")["clustervolume"].mean().sum() if data["size_cluster"] > 0 else 0.0
                assert data["bulk_volume"] >= data["bulk_clustervolume_cumulated"]
                data["bulk_density_free"] = df_nonsurface_free / (data["bulk_volume"] - data["bulk_clustervolume_cumulated"])
                data["surface_affinity_translation"] = df["surface_affinity_translation"].mean()
                data["surface_affinity_rotation"] = df["surface_affinity_rotation"].mean()

            sysdata.append(data)


        if clstr:
            for cluster_id, cluster in df.groupby("cluster"):
                if cluster["clustersize"].unique()[0] < clstr_min_size:
                    continue
                
                m = hashlib.md5()
                m.update(path_str_for_hash.encode('utf-8'))
                m.update(str(meta["time"]).encode('utf-8'))
                m.update(str(cluster_id).encode('utf-8'))
                identifier = m.hexdigest()
                
                data = {}
                data["id"] = identifier
                data["sys_id"] = sys_identifier
                data["simulation"] = path_id
                
                data["time"] = meta["time"]
                data["gamma"] = df["gamma"].mean()
                data["sigma"] = df["sigma"].mean() 
                data["kappa"] = df["kappa"].mean()
                data["epsilon"] = df["epsilon"].mean()
                data["x"] = meta["box.x"]
                data["y"] = meta["box.y"]
                data["z"] = meta["box.z"]
                data["volume"] = data["x"]*data["y"]*data["z"]
                data["temperature"] = np.round(meta["temperature"], 3)
                data["cluster_min_size"] = clstr_min_size
                
                data["x_mean"] = cluster["shiftx"].mean()
                data["y_mean"] = cluster["shifty"].mean()
                data["z_mean"] = cluster["shiftz"].mean()
                
                data["epot"] = cluster["epot"].mean()
                data["surfacepot"] = cluster["surfacepot"].mean()
                data["externalpot"] = cluster["externalpot"].mean()
                
                data["density"] = df.index.size / (meta["box.x"]*meta["box.y"]*meta["box.z"])
                data["density_free"] = float(df_free.index.size )/((meta["box.x"]*meta["box.y"]*meta["box.z"]) - (df_clstr.groupby("cluster")["clustervolume"].mean().sum() if df_clstr.index.size > 0 else 0.0))
                
                data["size"] = cluster.index.size

                data["chi"] = cluster["chi"].mean()
                data["order_simple"] = cluster["order"].mean()
                data["neighbors"] = cluster["neighbors"].mean()
                data["volume"] = cluster["clustervolume"].mean()

                clstrdata.append(data)

    return_dict = {}
    if sys:
        return_dict["sys"] = pd.DataFrame(sysdata)#.set_index('id')
    if clstr:
        return_dict["clstr"] = pd.DataFrame(clstrdata)#.set_index('id')
    
    return return_dict




def _insert_or_replace(df, cur, table_name, alter=True):
    type_map = {
        "int" : "INTEGER",
        "int32" : "INTEGER",
        "int64" : "INTEGER",
        "uint" : "INTEGER",
        "uint32" : "INTEGER",
        "uint64" : "INTEGER",
        "float" : "REAL",
        "float32" : "REAL",
        "float64" : "REAL",
        "object" : "TEXT",
        "bool" : "BOOL"
    }
    
    c_types = [type_map[t] for t in df.dtypes.values.astype(str)]
    c_names = ', '.join('{}'.format(n) for n in df.columns)
    c_names_at_creation = ', '.join('{} {}'.format(n,t) for n,t in zip(df.columns,c_types)).replace("id TEXT", "id TEXT PRIMARY KEY")
    c_names_at_creation = c_names_at_creation.replace("sys_id TEXT PRIMARY KEY", "sys_id TEXT")

    if not "system" in table_name:
        cur.execute(f"CREATE TABLE IF NOT EXISTS {table_name} ({c_names_at_creation}, FOREIGN KEY (sys_id) REFERENCES system_states (id))")
    else: 
        cur.execute(f"CREATE TABLE IF NOT EXISTS {table_name} ({c_names_at_creation})")

    if alter:
        db_columns = [i[1] for i in cur.execute(f"PRAGMA table_info({table_name})")]
        for df_column in df.columns:
            if df_column not in db_columns:
                cur.execute(f'ALTER TABLE {table_name} ADD COLUMN {df_column} {type_map[str(df[df_column].dtype)]} DEFAULT NULL')

    # df.loc[:, df.dtypes == object] = "'" + df.loc[:, df.dtypes == object]  + "'"

    for _, row in df.iterrows():
        values = ', '.join('{}'.format(k) for k in row.values)
        # cur.execute(f"INSERT OR REPLACE INTO {table_name} ({c_names}) VALUES ({values})".replace("nan", "NULL"))
        cur.executemany(f"INSERT OR REPLACE INTO {table_name} ({c_names}) VALUES ({', '.join('?' for _ in df.columns)})".replace("nan", "NULL"), df.values)



def dataToSQL(data, con):
    data["sys"].to_sql(name="system_states", con=con, if_exists="append")
    data["clstr"].to_sql(name="cluster_states", con=con, if_exists="append")
    # if "sys" in data.keys():
    #     try:
    #         if data["sys"].index.size > 0:
    #             cursor = con.cursor()
    #             try:
    #                 _insert_or_replace(data["sys"], cursor, "system_states")
    #             finally:
    #                 cursor.close()
    #     except Exception as e:
    #         print(f"Failed to write data.sys to sql connection: ", e)

    # if "clstr" in data.keys():
    #     try:
    #         if data["clstr"].index.size > 0:
    #             cursor = con.cursor()
    #             try:
    #                 _insert_or_replace(data["clstr"], cursor, "cluster_states")
    #             finally:
    #                 cursor.close() 
    #     except Exception as e:
    #         print(f"Failed to write data.clstr to sql connection: ", e)



def dataBufferToSQL(data_buffer, con):
    sys_df = pd.concat([d["sys"] for d in data_buffer if "sys" in d.keys()])
    print("writing", sys_df.index.size, "system states")
    sys_df.to_sql(name="system_states", con=con, if_exists="append")

    clstr_df = pd.concat([d["clstr"] for d in data_buffer if "clstr" in d.keys()])
    print("writing", clstr_df.index.size, "cluster states")
    clstr_df.to_sql(name="cluster_states", con=con, if_exists="append")

    # try:
    #     if sys_df.index.size > 0:
    #         cursor = con.cursor()
    #         try:
    #             _insert_or_replace(sys_df, cursor, "system_states")
    #         finally:
    #             cursor.close()
    # except Exception as e:
    #     print(f"Failed to write data.sys to sql connection: ", e)

    # try:
    #     if clstr_df.index.size > 0:
    #         cursor = con.cursor()
    #         try:
    #             _insert_or_replace(clstr_df, cursor, "cluster_states")
    #         finally:
    #             cursor.close() 
    # except Exception as e:
    #     print(f"Failed to write data.clstr to sql connection: ", e)



def calcProcessTime(starttime, cur_iter, max_iter):

    telapsed = time.time() - starttime
    testimated = (telapsed/cur_iter)*(max_iter)

    finishtime = starttime + testimated
    finishtime = dt.datetime.fromtimestamp(finishtime).strftime("%H:%M:%S")  # in time

    lefttime = testimated-telapsed  # in seconds

    return (int(telapsed), int(lefttime), finishtime)




def progress(status, remaining, total):
    print(f'Copied {total-remaining} of {total} pages...')



def gatherStates(
    basedir = Path("."),
    dbpath = "states.db",
    filenames = ["data.h5"],
    threads = 4,
    sys = True,
    clstr = True,
    key_prefix = "/time",
    clstr_min_size = 30,
    timeout = 30,
    buffer_size = 10
):
    # if not isinstance(basedir, type(Path("."))):
    if not isinstance(basedir, Path):
        basedir = Path(basedir)

    if not isinstance(dbpath, Path):
        dbpath = Path(dbpath)

    filepaths = []
    for fn in filenames:
        for fp in basedir.rglob(fn):
            filepaths.append(fp)
            # break

    NUM_FILES_FOUND = len(filepaths)
    print(f"found {NUM_FILES_FOUND} files")

    file_counter = 0
    
    print("database path:", dbpath)
    with sqlite3.connect(':memory:') as con:

        try:
            with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
                future_to_fp = { executor.submit(analyzeStates, path=fp, threads=1, key_prefix=key_prefix, sys=sys, clstr=clstr, clstr_min_size=clstr_min_size):fp for fp in filepaths}
                
                data_buffer = []

                START_TIME = time.time()
                for future in concurrent.futures.as_completed(future_to_fp):
                    file_counter += 1

                    fp = future_to_fp[future]
                    print()
                    print("[LOADED ] ", fp)
                    try:
                        # data = future.result(timeout=timeout)
                        data_buffer.append(future.result())
                    except Exception as exc:
                        print('%r generated an exception: %s' % (fp, exc))
                    # else:
                    #     dataToSQL(data, con)
                    #     print("[WRITTEN] ", fp)
                    #     runtime, left, eta = calcProcessTime(START_TIME, file_counter, NUM_FILES_FOUND)
                    #     print("[RUNTIME] ", f"time elapsed: {runtime:>6.0f} s, time left: {left:>6.0f} s, estimated finish time: {eta}")

                    if len(data_buffer) >= buffer_size:
                        print("[BUFFER ] ", "full")
                        dataBufferToSQL(data_buffer, con)
                        print("[BUFFER ] ", "written to", dbpath)
                        runtime, left, eta = calcProcessTime(START_TIME, file_counter, NUM_FILES_FOUND)
                        print("[RUNTIME] ", f"time elapsed: {runtime:>6.0f} s, time left: {left:>6.0f} s, estimated finish time: {eta}")
                        data_buffer = []
                        # gc.collect()
                
                if len(data_buffer):
                    print("[BUFFER ] ", "writing..")
                    dataBufferToSQL(data_buffer, con)
                    print("[BUFFER ] ", "written to", dbpath)
                    runtime, left, eta = calcProcessTime(START_TIME, file_counter, NUM_FILES_FOUND)
                    print("[RUNTIME] ", f"time elapsed: {runtime:>6.0f} s, time left: {left:>6.0f} s, estimated finish time: {eta}")
                    data_buffer = []
                    
        except TimeoutError as te:
            print(te)
        
        filedb = sqlite3.connect(dbpath)
        with filedb:
            con.backup(filedb)
        filedb.close()





def readStates(dbpath, sql):
    with sqlite3.connect(str(dbpath)) as con:
        return pd.read_sql_query(con=con, sql=sql)