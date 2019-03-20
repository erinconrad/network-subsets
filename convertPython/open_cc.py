def open_cc_f(file,whichsec):
    # import numpy
    import numpy as np

    # load the file
    data = np.load(file)

    # get list
    a = data[:,whichsec].tolist()
    
    return a


