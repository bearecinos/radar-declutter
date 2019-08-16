'''Allows one parameter of the model to be varied for multiple plots while
the others use the values defined in modelling.defaults (or passed to them).
Parameters to methods are generally the same as in radar.radargram(...)'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from scipy import signal
import multiprocessing as mp
import threading
from progress import progress
import h5py
import align
import modelling
from modelling.defaults import default
from modelling import waves, backscatter, directivity
import radar

env = modelling.parameters.env

# (wave, title) pairs
waveList = [(waves.GaussianDot(), "GaussianDot"), (waves.Ricker(), "Ricker"),
            (waves.CosExp(), "CosExp"), (waves.Sym(), "Sym")]

direcList = [(directivity.constant, "None"), (directivity.broad, "Broad"),
             (directivity.idl, "IDL"), (directivity.lobes, "Lobes")]

intensityList = [(backscatter.lambertian, "Diffuse"),
                 (backscatter.Min2, "Minnaert k=2"),
                 (backscatter.rayModel, "Ray tracing n=1"),
                 (backscatter.specular6, "cos^6")]


def _getFiles(name):
    try:
        files = os.listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x: int(x[5:-5]))  # assumes 'pointX.hdf5'
    return files


def compareWaves(name, adjusted=False, intensityModel=None, directional=None,
                 clip=0, rFactor=0.0, parallel=True):
    if intensityModel is None:
        intensityModel = default.getIntensity()
    if directional is None:
        directional = default.getDirectivity()

    print env
    plt.rcParams['axes.formatter.limits'] = [-4, 4]  # use standard form
    plt.figure(figsize=env.figsize)
    files = _getFiles(name)

    output = np.full((len(waveList), len(files), env.getSteps()), 0, float)

    cells = int(np.ceil(np.sqrt(len(waveList))))*110
    for i in range(len(waveList)):
        plt.subplot(cells+i+1)
        wave, title = waveList[i]
        output[i] = _radargram(name, files, adjusted, wave, intensityModel,
                               title, directional, clip, rFactor, parallel)
    plt.show()
    return output


def compareDirectivity(name, adjusted=False, intensityModel=None, wave=None,
                       clip=0, rFactor=0.0, parallel=True):
    if intensityModel is None:
        intensityModel = default.getIntensity()
    if wave is None:
        wave = default.getWave()

    print env
    plt.rcParams['axes.formatter.limits'] = [-4, 4]  # use standard form
    plt.figure(figsize=env.figsize)
    files = _getFiles(name)

    output = np.full((len(direcList), len(files), env.getSteps()), 0, float)

    cells = int(np.ceil(np.sqrt(len(direcList))))*110
    for i in range(len(direcList)):
        plt.subplot(cells+i+1)
        directional, title = direcList[i]
        output[i] = _radargram(name, files, adjusted, wave, intensityModel,
                               title, directional, clip, rFactor, parallel)
    plt.show()
    return output


def compareBackscatter(name, adjusted=False, wave=None, directional=None,
                       clip=0, rFactor=0.0, parallel=True):
    if directional is None:
        directional = default.getDirectivity()
    if wave is None:
        wave = default.getWave()

    print env
    plt.rcParams['axes.formatter.limits'] = [-4, 4]  # use standard form
    plt.figure(figsize=env.figsize)
    files = _getFiles(name)

    output = np.full((len(intensityList), len(files), env.getSteps()),
                     0, float)

    cells = int(np.ceil(np.sqrt(len(intensityList))))*110
    for i in range(len(intensityList)):
        plt.subplot(cells+i+1)
        intensityModel, title = intensityList[i]
        output[i] = _radargram(name, files, adjusted, wave, intensityModel,
                               title, directional, clip, rFactor, parallel)
    plt.show()
    return output


def _radargram(name, files, adjusted=False, wave=None, intensityModel=None,
               title=None, directional=None, clip=0, rFactor=0.0,
               parallel=True):

    returnData = np.full((len(files), env.getSteps()), 0, float)

    p = mp.Pool(mp.cpu_count())
    # data needed to run across several processors
    data = [(i, name+"/"+files[i], wave, intensityModel, directional, env,
             rFactor) for i in range(len(files))]
    try:
        if parallel:
            for i, ar in progress(p.imap_unordered(radar.worker, data),
                                  len(files)):
                returnData[i] = ar
        else:  # non-parallel option for full error reporting
            for j in progress(range(len(files))):
                i, ar = radar.worker(data[j])
                returnData[i] = ar
    except IOError as e:
        p.close()
        print "\nError reading hdf5 file :\n"+e.message
        return -1
    p.close()

    if adjusted:  # align first responses with top of radargram
        returnData = align.minAlign(returnData, env.getDx(), 200.0)

    ys = np.linspace(0, env.getMaxTime(), env.getSteps())

    if clip:  # clipping if set as non-zero
        returnData = np.clip(returnData, np.percentile(returnData, clip),
                             np.percentile(returnData, 100-clip))

    plt.ylim(env.getMaxTime(), 0)
    draw = np.swapaxes(returnData, 0, 1)

    plt.contourf(np.arange(len(files)), ys, draw, 100,
                 norm=radar.MidNorm(np.mean(draw)), cmap="Greys")
    if title is not None:
        plt.title(title)
    plt.colorbar()

    return returnData
