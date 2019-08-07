import numpy as np
# each model is a method of type float * float -> float
# should work with numpy arrays
# must use def rather than lambdas in order to pass with multiprocessing

def constant(theta,phi):
    return 1
def broad(theta,phi):
    return abs(np.sin(theta))**0.5
def idl(theta,phi):
    return np.sin(theta)**6
def lobes(theta,phi):
    return np.sin(theta)**2*np.sin(3*theta)**2
