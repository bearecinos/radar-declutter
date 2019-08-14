import numpy as np

# each model is a method of type float -> float
# should work with numpy arrays
# have to use def rather than lambdas in order to pass with multiprocessing

def lambertian(theta):
    return np.cos(theta)

def Minnaert(theta,k): 
    return np.power(np.cos(theta),2*k-1)
def Min2(theta):
    return Minnaert(theta,2)

def raySpecular(theta,n=1):
    return np.power(np.maximum(0,np.cos(theta*2)),n)
def rayModel(theta):
    return np.maximum(0,np.cos(theta*2))
def ray2Model(theta):
    return np.power(np.maximum(0,np.cos(theta*2)),2)

def specular2Model(theta):
    return np.power(np.cos(theta),2)
def specular6Model(theta):
    return np.power(np.cos(theta),6)
def IDLModel(theta):
    return np.power(np.cos(theta),8)
