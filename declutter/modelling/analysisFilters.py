'''Provides a collection of methods to mask points to consider as part of the
walls of a glacier. Performance can vary significantly with data so worth
testing. e.g. glacier appears very rough with some data, smooth with other.

Only interested in visible points so, depending on type of filter, both 1D or
2D arrays can be returned, and these are dealt with by compose method.'''
import numpy as np
from scipy import signal


def smooth(grid, threshold=0.95, sampleSteps=5, sampleArea=10):
    '''Where glacier is rougher than walls in data, this filters out
    sufficiently rough points. Note that the edges of the array will have
    poor accuracy.

    Parameters
    grid - 2D float arrray : The heightmap for the area.
    threshold - float : A value between 0 and 1. The correlation coefficient
        of the area sampled must be above this to be considered smooth enough.
    sampleSteps - int : How regularly to sample the smoothness of the array.
        This is a slow process so larger steps increase speed at the cost of
        accuracy.
    sampleArea - int : The side length of the square to consider the smoothness
        over.

    Returns
    mask - 2D bool array : True cells are ones which are sufficiently
        smooth.'''
    result = np.full((grid.shape[0]/sampleSteps, grid.shape[1]/sampleSteps),
                     0.0, float)
    ys, xs = np.indices((sampleArea, sampleArea))
    ys, xs = ys.reshape(-1), xs.reshape(-1)
    results = []
    for x in range(0, grid.shape[1]-sampleArea, sampleSteps):
        for y in range(0, grid.shape[0]-sampleArea, sampleSteps):
            region = grid[y:y+sampleArea, x:x+sampleArea].reshape(-1)
            if np.sum(np.isnan(region)) > 0:  # ignore invalid regions
                continue
            A = np.vstack([xs, ys, np.ones(len(region))]).T
            # least squares estimate
            residuals = np.linalg.lstsq(A, region)[1]
            # R^2 correlation coefficient
            results = 1.0 - residuals/np.sum(
                region-np.mean(region)**2)
    return _scaleUp(grid, results, sampleSteps, sampleArea) > threshold


def _scaleUp(grid, results, steps, area):
    '''Resamples the array to give a value for every cell based on the sampled
    area the point is the closest to the centre of. Points on the edges are
    furthest from the centre of a sample.'''
    height, width = results.shape
    y, x = np.indices(grid.shape)
    y = np.clip((y-area/2)/steps, 0, height-1)
    x = np.clip((x-area/2)/steps, 0, width-1)
    return results[y, x]


def steepness(slopeGrid, angle=45):
    '''Filters out flat surfaces which clearly aren't glacier wall.'''
    return slopeGrid > angle


def visibility(visible, threshold=1.0, sampleArea=5):
    '''This filters out cells which do not have all of the surrounding cells
    visible. This can help to ignore the glacier where the surface is rough.

    Parameters
    ----------
    visible - 2D bool array : The mask of visible cells
    threshold - float : The proportion of cells which must be visible, between
        0 and 1.
    sampleArea - int : The side length of the square to sum the number of
        visible cells in.

    Returns
    -------
    mask - 2D bool array : Mask of sufficiently visible points.
    '''
    return (signal.convolve(visible, np.full((5, 5), 1, int), "same") >
            threshold*sampleArea**2)


def height(grid, z, offset=20.0):
    return grid > z + offset


def glacierDirection(grid, corner, x, y, cellsize, cutOut=30):
    '''Uses a linear fit of the points within 50m altitude of the surface
    below the radar to estimate glacier direction.'''
    coords = [(x-corner[0])/cellsize, (y-corner[1])/cellsize]
    groundHeight = grid[int(coords[1]), int(coords[0])]
    ys, xs = np.indices(grid.shape)
    m = (grid > groundHeight - 50) & (grid < groundHeight + 50)
    glacierGrad = np.arctan(np.polyfit(xs[m], ys[m], deg=1)[0])
    # NaN points should be ignored (directly below radar)
    # inf values are correctly converted by arctan
    grad = (ys-coords[1]) / (xs-coords[0])

    # remove returning the angle
    return abs(np.arctan(grad)-glacierGrad) > cutOut * np.pi/180


def horizon(grid, visible, distance, z, cutOut=5):
    '''Masks out point which, from the radar, are at an angle greater than
    cutOut from the horizon.'''
    return abs(z-grid[visible])/distance < np.sin(cutOut*np.pi/180)


def minDist(distance, minimum=400):
    return distance >= minimum


def compose(filters, visible):
    """Return a mask which is true only for the points satisfying every mask
    given. Points must also be visible.
    Takes both 2D and 1D arrays, returns 1D mask to apply to visible points."""
    if len(filters) == 0:
        return visible
    result = np.full(visible.shape, True, bool)
    for f in filters:
        if len(f.shape) == 2:
            result &= f
    result = result[visible]
    for f in filters:
        if len(f.shape) == 1:
            result &= f
    return result
