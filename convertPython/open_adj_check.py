def open_ad_f(file):
    # import numpy
    import numpy as np
    # load the file
    data = np.load(file)
    keys = data.item().keys()
    vals = [];
    for k in keys:
        vals.append(data.item()[k])
    vals, b = zip(*vals)
    keys = zip(*sorted(zip(vals,keys)))[1]
    vals = sorted(vals)
    return vals, keys


