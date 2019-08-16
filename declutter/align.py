'''Used to shift each radargram slice so that the first response is at the
top of the plt. Where the first response is not detected within a given range,
the data is assumed to be incomplete and this point is interpolated from
neighbouring values.'''
import numpy as np


def which(x):  # 1D where simplification
    return np.where(x)[0]


def fillGaps(indices):
    """Given an array of values with NaNs in some positions, interpolate the
    valid points to replace the NaNs. NaNs at the end of the range are set to
    the first valid point closest to that end of the array."""
    nans = np.isnan(indices)  # mask of NaN values

    # interpolate the values for NaN indices given the values of
    # non-NaN values. Fill the NaN indices with these values.
    indices[nans] = np.interp(which(nans), which(~nans), indices[~nans])
    return indices.astype(int)


def minAlign(data, dx=1.875, cutoff=200.0):
    """Shifts the data for each point so the the first response for each
    point is at the start of the data.

    Parameters
    ----------
    data - 2D float array : The data to plot for each point.
    dx - float (optional) : The spatial distance between cells. Default 1.875m.
    cutoff - float (optional) : The distance at which to assume the surface
        directly below the radar was not seen, hence ignore the first
        response seen.
    """
    out = np.full_like(data, 0.0)
    height, width = data.shape
    indices = np.full(height, 0.0)
    for i in range(len(data)):
        w = np.where(data[i] != 0)[0]  # array of non-zero values
        if len(w) == 0 or w[0]*dx > cutoff:  # all 0 up to cutoff
            indices[i] = np.nan  # will be interpolated
        else:
            indices[i] = w[0]
    indices = fillGaps(indices)  # interpolate
    for i in range(len(data)):
        out[i, :width-indices[i]] = data[i, indices[i]:]  # shift in output
    return out
