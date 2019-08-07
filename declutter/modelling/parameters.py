import os
import numpy as np

# Just scales waves, all have own actual peak frequencies
wavelength = 100
freq = 3e8/wavelength # 3 MHz

class Env:
    maxDist = 3000.0
    maxTime = maxDist/1.5e8
    dt = 1.25e-8
    dx = dt*1.5e8
    steps = int(maxTime / dt)
    
env = Env()

def loadParameters():
    path = os.path.dirname(__file__)+"/config.npy"
    if not os.path.exists(path):
        return -1
    setups = {"maxDist":setMaxDist, "maxTime":setMaxTime, "dx":setSpaceStep,
              "dt":setTimeStep, "steps":setSteps}
    data = np.load(path,allow_pickle=True).item()
    print "Loading plot parameters from config file:"
    # Calling in certain orders changes some values back, hence cases
    if data["steps"] is None or (data["maxDist"] is None and data["maxTime"] is None):
        for key, val in data.iteritems():
            if val is not None:
                setups[key](val)
                print key+" : "+str(val)
    elif data["maxDist"] is not None:
        setMaxDist(data["maxDist"])
        setSpaceStep(data["maxDist"]/data["steps"])
        print "maxDist : "+str(data["maxDist"])
        print "steps : "+str(data["steps"])
    else:
        setMaxTime(data["maxTime"])
        setSpaceStep(data["maxTime"]/data["steps"])
        print "maxTime : "+str(data["maxTime"])
        print "steps : "+str(data["steps"])
    return 0

def storeParameters(env = env):
    # Rest found from maxDist and dt on loading
    data = {"maxDist":None,"maxTime":env.maxTime,
            "dx":None, "dt":env.dt, "steps":None}
    path = os.path.dirname(__file__)+"/config.npy"
    np.save(path,data)
    return 0

# Distance is for speed in air rather than ice. i.e. 3e8, not 1.67e8
def setTimeStep(dt = 1.25e-8):
    """Sets the time between sample points in radargrams/wiggle plots.
    This is in terms of the two-way path i.e. difference in recieved time.
    Default is 1.25e-8s"""
    env.dt = dt
    env.dx = dt*1.5e8
    _setSteps()
    
def setSpaceStep(dx = 1.875):
    """Sets the distance between sample points in the radargram/wiggle plots.
    This is in terms of the one-way path. Default is 1.875m"""
    env.dx = float(dx)
    env.dt = dx/1.5e8
    _setSteps()
    
def setMaxDist(d = 3000.0):
    """Sets the maximum range to consider a surface creating a response from.
    Default is 3000m."""
    env.maxDist = float(d)
    env.maxTime = d/1.5e8
    _setSteps()
    
def setMaxTime(t = 2e-5):
    """Sets the maximum duration to show the radargram/wiggle plot over.
    Default is 2e-5s."""
    env.maxTime = t
    env.maxDist = t*1.5e8
    _setSteps()
    
def setSteps(n = 1600):
    """Sets the number of samples in the radargram/wiggle plot. This
    changes the maximum range of the plot. Default is 1600."""
    env.steps = int(n)
    env.maxDist = env.dx*n
    env.maxTime = env.dt*n
    
def _setSteps():
    env.steps = int(env.maxTime / env.dt)

figsize = (12,6.8)
def setFigSize(x,y):
    """Sets the size of any plots. Default is (12,6.8)"""
    global figsize
    figsize = (x,y)
