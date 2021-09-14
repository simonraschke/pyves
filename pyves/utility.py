import pandas as pd
import numpy as np
import json



def h5store(filename, key, df, **kwargs):
    from pathlib import Path
    store = pd.HDFStore(filename, mode="a", complevel=1)
    store.put(key, df, format="table")
    store.get_storer(key).attrs.metadata = kwargs
    store.close()



def _h5load_inner(store, key):
    data = store[key]
    try:
        metadata = store.get_storer(key).attrs.metadata
    except Exception as e:
        raise Warning("No metadata")
        return data, None
    return data, metadata



def h5load(filename, key):
    with pd.HDFStore(filename) as store:
        data, metadata = _h5load_inner(store, key)
    return data, metadata



def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))



def generateSystem(datapath, parameterpath, time):
    df, metadata = h5load("test/data.h5", "time1500")
    # _prms = pyves.parametersToDict("parametrs.json")
    with open("test/parameters.json", 'r') as prms_file:
        _prms = json.loads(prms_file.read())
    sysfromdata = _makeSystem(df, _prms, metadata)
    return sysfromdata



def _makeSystem(
        df : pd.DataFrame, 
        prms : dict,
        metadata : dict
    ):
    assert("x" in df.columns)
    assert("y" in df.columns)
    assert("z" in df.columns)
    assert("ux" in df.columns)
    assert("uy" in df.columns)
    assert("uz" in df.columns)
    assert isinstance(prms, dict)
    assert isinstance(metadata, dict)

    import _pyves
    system = _pyves.System()
    assert isinstance(system, _pyves.System)
    system.threads = 1 #prms["hardware"].get("threads", 1)
    system.box.x = metadata["box.x"]
    system.box.y = metadata["box.y"]
    system.box.z = metadata["box.z"]
    system.temperature = metadata["temperature"]
    system.interaction_cutoff = metadata["interaction.cutoff"]
    system.cell_update_interval = prms["control"]["cell_update_interval"]
    system.neighbor_update_interval = prms["control"]["neighbor_update_interval"]
    system.neighbor_cutoff = prms["control"]["neighbor_cutoff"]
    system.interaction_surface = bool(metadata["interaction.surface"])
    system.interaction_surface_width = metadata["interaction.surface_width"]

    for _, row in df.iterrows():
        system.particles.append(_pyves.Particle([row["x"], row["y"], row["z"]], [row["ux"], row["uy"], row["uz"]], 
            sigma=row["sigma"], kappa=row["kappa"], eps=row["epsilon"], gamma=row["gamma"], name=row["name"]))
        system.particles[-1].surface_affinity_translation = row["surface_affinity_translation"]
        system.particles[-1].surface_affinity_rotation = row["surface_affinity_rotation"]
    
    # Setup cells
    box_dims = np.array([system.box.x, system.box.y, system.box.z])
    cells_per_dim = np.array(box_dims/metadata["cell_min_size"]).astype(int)
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
   
    return system