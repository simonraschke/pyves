import pandas as pd



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