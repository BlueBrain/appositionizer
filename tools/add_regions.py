"""Small utility to edit circuit files
"""
import sys

import h5py as h5
import numpy as np

filename, = sys.argv[1:]

cnt = 0

with h5.File(filename, 'r+') as f:
    layers = f['/cells/properties/layer']
    data = np.zeros_like(layers)
    for i in range(len(data)):
        if i < 8:
            data[i] = 5
        elif i % 3 == 0 and i % 5 == 0:
            if cnt < 50:
                data[i] = 4
                cnt += 1
            else:
                data[i] = 3
        elif i % 3 == 0:
            data[i] = 1
        elif i % 5 == 0:
            data[i] = 2
    dt = h5.special_dtype(vlen=bytes)
    try:
        del f['/cells/properties/region']
        del f['/library/region']
    except:
        pass
    f.create_dataset('/cells/properties/region', data=data, dtype='<i4')
    f.create_dataset('/library/region', data=np.array([b'L1', b'L2', b'L3', b'L23', b'L123', b'L456']), dtype=dt)
