import numpy as np

def fillGaps(indices):
    nans = np.isnan(indices)
    which = lambda x : np.where(x)[0]
    indices[nans] = np.interp(which(nans),which(~nans),indices[~nans])
    return indices.astype(int)

def minAlign(data, dx=1.875, cutoff=200.0):
    out = np.full_like(data,0.0)
    height,width = data.shape
    indices = np.full(height,0.0)
    for i in range(len(data)):
        w = np.where(data[i] != 0)[0]
        if len(w) == 0 or w[0] > cutoff:
            indices[i] = np.nan
        else:
            indices[i] = w[0]
    indices = fillGaps(indices)
    for i in range(len(data)):
        out[i,:width-indices[i]] = data[i,indices[i]:]
    return out
