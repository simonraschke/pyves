
import json
import collections



# https://stackoverflow.com/a/30655448
def deep_update(source, overrides):
    """
    Update a nested dictionary or similar mapping.
    Modify ``source`` in place.
    """
    for key, value in overrides.items():
        if isinstance(value, collections.Mapping) and value:
            returned = deep_update(source.get(key, {}), value)
            source[key] = returned
        else:
            source[key] = overrides[key]
    return source



_particle_string = \
"""
{
    "RAND" :
    {
        "number" : 175,
        "dist" : "random",
        "sigma" : 1.0,
        "kappa" : 1.0,
        "epsilon" : 1.0,
        "gamma" : 0.1914626,
        "bound_translation" : null,
        "bound_rotation" : null,
        "surface_affinity_translation" : 0,
        "surface_affinity_rotation" : 0,
        "self_affinity" : 1.01,
        "other_affinity" : 0.99
    }
}
"""
default_particle = json.loads(_particle_string)



_system_str = \
"""
{
    "temperature" : 0.20,
    "global_exchange_ratio" : 0.05,
    "global_exchange_epot_theshold" : -1.0,
    "box" : 
    {
        "x" : 20,
        "y" : 20,
        "z" : 20
    },
    "translation" :
    {
        "min" : 0.2,
        "max" : 0.2,
        "target" : 0.3,
        "interval" : 1000
    },
    "rotation" :
    {
        "min" : 0.2,
        "max" : 0.2,
        "target" : 0.3,
        "interval" : 1000
    },
    "interaction" : 
    {
        "cutoff" : 3.0,
        "z_surface" : false,
        "z_surface_width" : 1.0
    },
    "particles" : {
    }
}
"""
default_system = json.loads(_system_str)
deep_update(default_system, {"particles":default_particle})



_input_string = \
"""
{
    "hardware" : 
    {
        "threads" : 4
    },
    "system" : {
    },
    "control" : 
    {
        "time_delta" : 1000,
        "time_max" : 1500,
        "cell_min_size" : 3.5,
        "cell_update_interval" : 3,
        "neighbor_update_interval" : 3,
        "neighbor_cutoff" : 3.5
    },  
    "output" :
    {
        "dir" : "test",
        "filename" : "data.h5",
        "mode" : "overwrite",
        "direct_analysis" : false
    },
    "input" :
    {
        "dir" : ".",
        "filename" : "data.h5",
        "key" : "HEAD"
    },
    "analysis":
    {
        "volume_max_grid" : 30,
        "volume_points_per_sigma" : 3,
        "DBSCAN_eps" : 1.2,
        "neighbor_cutoff" : 1.2
    }
}
"""
default_prms = json.loads(_input_string)
deep_update(default_prms, {"system":default_system})