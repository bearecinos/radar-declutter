import os
import numpy as np

# Just scales waves, all have own actual peak frequencies
wavelength = 100
freq = 3e8/wavelength # 3 MHz

class Env:
    dt = 1.25e-8
    steps = 1600
    def __str__(self):
        return "Settings: maxTime = "+str(self.dt*self.steps)+", dt = "+str(self.dt)+", steps = "+str(self.steps)
    def getDt(self):
        return self.dt
    def getDx(self):
        return self.dt*1.5e8
    def getSteps(self):
        return self.steps
    def getMaxTime(self):
        return self.steps * self.dt
    def getMaxDist(self):
        return self.steps * self.dt * 1.5e8

# unclear how to in python, but should only ever have 1 instance created
# e.g. singleton design pattern
env = Env() # parameters loaded at end of module

def loadParameters():
    global env    
    path = os.path.dirname(__file__)+"/config.npy"
    if not os.path.exists(path):
        return -1
    data = np.load(path,allow_pickle=True).item()

    # will always overwrite anything previously set
    setTimeStep(data["dt"])
    setSteps(data["steps"])
    return 0

def storeParameters(env = env):
    # Rest found from maxDist and dt on loading
    data = {"dt":env.dt, "steps":env.steps}
    path = os.path.dirname(__file__)+"/config.npy"
    np.save(path,data)
    return 0

# Distance is for speed in air rather than ice. i.e. 3e8, not 1.67e8
def setTimeStep(dt = 1.25e-8):
    """Sets the time between sample points in radargrams/wiggle plots.
    This is in terms of the two-way path i.e. difference in recieved time.
    Default is 1.25e-8s. Maximum range unchanged."""
    env.steps = int(env.steps * (env.dt / dt) + 0.5) # near rounding
    env.dt = dt
    
def setSpaceStep(dx = 1.875):
    """Sets the distance between sample points in the radargram/wiggle plots.
    This is in terms of the one-way path. Default is 1.875m. Maximum range unchanged."""
    setTimeStep(dx/1.5e8)
   
def setMaxTime(t = 2e-5):
    """Sets the maximum duration to show the radargram/wiggle plot over.
    Default is 2e-5s. Time granularity unchanged."""
    env.steps = int(t/env.dt + 0.5) # near rounding
    
def setMaxDist(d = 3000.0):
    """Sets the maximum range to consider a surface creating a response from.
    Default is 3000m. Time granularity unchanged."""
    setMaxTime(d/1.5e8)
    
def setSteps(n = 1600):
    """Sets the number of samples in the radargram/wiggle plot. This
    changes the maximum range of the plot. Default is 1600."""
    env.steps = int(n + 0.5) # near rounding
 
loadParameters()

figsize = (12,6.8)
def setFigSize(x,y):
    """Sets the size of any plots. Default is (12,6.8)"""
    global figsize
    figsize = (x,y)
    
