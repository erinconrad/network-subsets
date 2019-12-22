# import scipy
import scipy.io as io
# import numpy
import numpy as np
import sys

files = [('alphatheta', 'all_adj_alphatheta.npy'),
             ('beta', 'all_adj_beta.npy'),
             ('broadband', 'all_adj_broadband_CC.npy'),
             ('highgamma', 'all_adj_highgamma.npy'),
             ('lowgamma', 'all_adj_lowgamma.npy'),
             ('veryhigh', 'all_adj_veryhigh.npy'),
             ('labels', 'labels.npy')]
output = 'all.mat'

def save_as_mat_type1(dir):
    out = {}
    for name, path in files:
        d = np.load(dir + '/' + path)
        if name == 'labels':
            keys = d.item().keys()
            vals = [];
            for k in keys:
                vals.append(d.item()[k])
            vals, b = zip(*vals)
            out['label_vals'] = vals
            out['label_keys'] = keys
            out['label_b'] = b
        else:
            out[name] = d
    io.savemat(dir + '/' + output, out)

if __name__ == '__main__':
    #print('here are the arguments I got:')
    #for i, arg in enumerate(sys.argv):
    #    print('#{}: {}'.format(i, arg))
    save_as_mat_type1(sys.argv[1])
