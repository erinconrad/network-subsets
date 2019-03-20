def get_shape_f(file):
    # import numpy
    import numpy as np

    # load the file
    data = np.load(file)
    
    return data.shape
