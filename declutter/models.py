"""Uses a 'maps.hdf5' file and directory of .hdf5 files for points to model
the radargram/wiggle plots seen. Also contains a range of methods to allow
the model to be altered, such as the backscatter model or the wave to convolve
with the result."""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from os import listdir
from scipy import signal
import multiprocessing as mp
from progress import progress
import h5py


__all__ = ["setFigSize", "setTimeStep", "setSpaceStep", "setMaxDist", "setMaxTime",
           "setSteps", "lambertian", "Minnaert", "Min2", "raySpecular", "rayModel",
           "ray2Model", "specular2Model", "specular4Model", "IDLreflection",
           "GaussianDot", "Ricker", "Constant", "IDLWave", "RC", "CosExp", "Sym",
           "direcBroad", "direcNone", "direcIDL", "direcLobes", "loadArrays",
           "processSlice", "models", "titles", "compare", "wiggle", "manyWiggle",
           "showWave", "showDirectionality"]

figsize = (12,6.8)
def setFigSize(x,y):
    """Sets the size of any plots. Default is (12,6.8)"""
    figsize = (x,y)

# Distance is for speed in air rather than ice. i.e. 3e8, not 1.67e8
def setTimeStep(dt = 1.25e-8):
    """Sets the time between sample points in radargrams/wiggle plots.
    This is in terms of the two-way path i.e. difference in recieved time.
    Default is 1.25e-8s"""
    global env
    env.dt = dt
    env.dx = dt*1.5e8
    _setSteps()
    
def setSpaceStep(dx = 1.875):
    """Sets the distance between sample points in the radargram/wiggle plots.
    This is in terms of the one-way path. Default is 1.875m"""
    global env
    env.dx = float(dx)
    env.dt = dx/1.5e8
    _setSteps()
    
def setMaxDist(d = 3000.0):
    """Sets the maximum range to consider a surface creating a response from.
    Default is 3000m."""
    global env
    env.maxDist = float(d)
    env.maxTime = d/1.5e8
    _setSteps()
    
def setMaxTime(t = 2e-5):
    """Sets the maximum duration to show the radargram/wiggle plot over.
    Default is 2e-5s."""
    global env
    env.maxTime = t
    env.maxDist = t*1.5e8
    _setSteps()
    
def setSteps(n = 1600):
    """Sets the number of samples in the radargram/wiggle plot. This
    changes the maximum range of the plot. Default is 1600."""
    global env
    env.steps = n
    env.maxDist = dx*n
    env.maxTime = dt*n
    
def _setSteps():
    global env
    env.steps = int(env.maxTime / env.dt)

class Env:
    maxDist = 3000.0
    maxTime = maxDist/1.5e8
    dt = 1.25e-8
    dx = dt*1.5e8
    steps = int(maxTime / dt)
    
env = Env()
    

def _time(distance):
    """Get the time for a wave to return to the radar after reflecting off a
    point 'distance' away."""
    return distance*2.0/3e8

# surface models
def lambertian(theta):
    return np.cos(theta)
def Minnaert(theta,k): 
    return np.power(np.cos(theta),2*k-1)
def Min2(theta):
    return Minnaert(theta,2) # as can't pass lambda to processor pool
def raySpecular(theta,n=1):
    return np.power(np.maximum(0,np.cos(theta*2)),n)
def rayModel(theta):
    return np.maximum(0,np.cos(theta*2))
def ray2Model(theta):
    return np.power(np.maximum(0,np.cos(theta*2)),2)
def specular2Model(theta):
    return np.power(np.cos(theta),2)
def specular4Model(theta):
    return np.power(np.cos(theta),4)
def IDLreflection(theta):
    return np.power(np.cos(theta),6)

##########
_wavelength = 100
_freq = 3e8/_wavelength

def bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return signal.lfilter(b, a, data)

class MidNorm(colors.Normalize):
    """ Normalize colorbar so that diverging bars work from chosen value"""
    def __init__(self,midpoint,vmin=None,vmax=None,clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self,vmin,vmax,clip)
    def __call__(self,value,clip=None):
        x,y = [self.vmin,self.midpoint,self.vmax],[0,0.5,1]
        return np.ma.masked_array(np.interp(value,x,y),np.isnan(value))

# wave models - some rely on environment variables
class GaussianDot:
    delt = 2.0/(3.0*_freq)
    align = 0.0
    def amplitude(self,t):
        return (-np.exp(0.5)*2*np.pi*_freq*(t-self.delt)*
                np.exp(-2*np.pi**2*_freq**2*(t-self.delt)**2))
# looks poor for one dataset but good for other?
class Ricker:
    delt = 2.0/(3.0*_freq)
    align = 0.0
    def amplitude(self,t):
        t = t - self.delt
        return (1-2*np.pi*np.pi*t*t*_freq*_freq)*np.exp(
            -np.pi*np.pi*_freq*_freq*t*t)

class Constant:
    delt = 0.0
    align = 0.0
    def amplitude(self,t):
        return t < env.dt 
class IDLWave:
    # exp( - (findgen(echo_size)-float(i))^2/50. )
    # timestep was 240us (2.4e-4 seconds echo length) / 2048 granularity
    # = 1.172e-7 = 117ns - much longer than our sampling rate
    # scaling = 1.172e-7
    align = env.maxTime / 2.0
    delt = align
    def __init__(self,c=50.0):
        self.c = float(c)
    def amplitude(self,t):
        t = (t-self.delt)/env.dt
        return np.exp(-t**2/self.c)
class RC:
    align = 0.0
    def __init__(self,c=5.0):
        self.c = float(c)
    def amplitude(self,t):
        result = np.full_like(t,0.0) # get around vectorizing
        result[t>0] = np.exp(-t[t>0]/(self.c*1e-6))
        return result
    
class CosExp:
    delt = 1.1/_freq
    align = 0.0
    def amplitude(self,t):
        p = 2*np.pi*_freq*(t-self.delt)
        return np.cos(1.5*p)*np.exp(-p*p/16.0)
class Sym:
    delt = 0.5/_freq
    align = 0.0
    def amplitude(self,t):
        p = 2*np.pi*_freq*(t-self.delt)
        return (160*p-56*p*p*p+3*p**5)*np.exp(-p*p/4)/(1.343*64.0)


def direcNone(theta,phi):
    return 1
def direcBroad(theta,phi):
    return abs(np.sin(theta))**0.5
def direcIDL(theta,phi):
    return np.sin(theta)**6
def direcLobes(theta,phi):
    return np.sin(theta)**2*np.sin(3*theta)**2

def loadArrays(filename):
    """Loads the data for a point. Will raise an IOError if the file
    does not exist or has the wrong form."""
    with h5py.File(filename,"r") as f:
        distance = f["distance"][()]
        angle = f["incidence"][()]
        theta,phi = None, None
        if "antennaTheta" in f:
            theta = f["antennaTheta"][()]
            phi = f["antennaPhi"][()]
    return distance,angle,theta,phi

# Replace GaussianDot with RC
def processSlice(distance,angle,theta,phi,intensityModel=raySpecular,
                 wave=GaussianDot(),rFactor=0,directional=direcNone):
    """Models the response seen at a single point.

    Parameters
    distance - float array : the distance in metres to every visible point.
    angle - float array : the surface incidence angle (in radians) of every
        visible point.
    theta, phi - float arrays : the direction to each visible point in spherical
        coordinates, with the ends of the antenna being the poles. These can
        be None, in which case they are not considered.
    intensityModel - function (optional) : a function of incidence angle for the
        backscatter from a surface. raySpecular by default i.e. cos(2 theta)
    wave - class instance (optional) : a model of the wave to convolve with the
        response. GaussianDot() by default.
    rFactor - float (optional) : How intensity should fall with distance. 0 by default.
    directional - function  (optional) : A function of theta and phi for the
        directivity of the antenna.
        
    Returns
    convol - float array : A time series of the modelled response.
    """
    m = (distance < env.maxDist)
    t = (_time(distance[m])/env.dt).astype(int)
    
    intensity = intensityModel(angle[m]) / np.power(distance[m],rFactor)
    
    # directionality 
    # input function to use
    if theta is not None:
        intensity *= directional(theta[m],phi[m])
    
    sample = np.array([np.sum(intensity[t==i]) for i in range(env.steps)])
    #sample = np.array([np.sum(t==i) for i in range(_steps)])
    #sample = np.array([np.sum(intensity[t==i])/(1+np.sum(t==i)) for i in range(_steps)])
    
    # convolution
    w = wave.amplitude(np.linspace(0.0,env.maxTime,env.steps))
    idx = int(wave.align/env.dt)
    convol = np.convolve(sample,w)[idx:idx+env.steps]
    
    # Filtering
    #convol = bandpass_filter(convol,0.05,0.4,1.0)

    return convol


models = [rayModel]#,ray2Model,specular2Model,specular4Model]
titles = ["Ray tracing n=1"]#,"Ray tracing n=2","sin(theta)^2","sin(theta)^4"]

def compare(name,adjusted=False,wave=GaussianDot(),save=None,models=models,
            titles=titles,directional=direcNone,display=True):
    '''Plots the radargram for the points in the given directory once for each
    model given in the models list. Adjusted aligns the y-axis by elevation
    rather than timing. The wave used by the model can also be changed and setting
    "save" to a string means the plot will be saved to that location.

    Parameters
    name - string : directory to look in for point data. Must contain no other
                    files or directories.
    adjusted - bool : indicates that the responses from each radar point should
                      be aligned vertically. i.e. If the radar took two samples
                      at the same coordinates but different elevations, the
                      response from the surface directly below would be at the
                      same point on the plot for both samples.
    wave - Wave instance : the wave to convolve over the response.
    save - string (optional) : the name of the file to save the radargram in.
    models - a list of functions of incidence angle to use for modelling the response.
    titles - a list of names to display above each plot, corresponding to the model
        in the same index in the models list.

    Returns
    returnData - 2D float array : The generated radargram data.
    Returns -1 if the method fails.
    '''
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=figsize)
    try:
        files = listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes 'pointX.hdf5'
    heights = []
    returnData = np.full((len(models),len(files),env.steps),0,float) # 3D - many plots
    
    p = mp.Pool(mp.cpu_count())
    data = [(i,name+"/"+files[i],wave,models,directional,env) for i in range(len(files))]

    try:
        for i, h, ars in progress(p.imap_unordered(worker,data),len(files)):
            returnData[:,i] = ars
            heights.append(h)
    except IOError as e:
        p.close()
        print "\nError reading hdf5 file :\n"+e.message
        return -1
    p.close()

    highest,lowest = 0,0
    if adjusted:
        highest = max(heights)
        lowest = min(heights)
        draw = np.full((len(models),len(files),env.steps + int((highest-lowest)/env.dx)),0)
        for i in range(len(files)):
            start = int((highest-heights[i])/env.dx)
            draw[:,i,start:start+env.steps] = returnData[:,i]
        returnData = draw

    cells = int(np.ceil(np.sqrt(len(models))))*110
    ys = np.linspace(0,(env.maxDist+highest-lowest)*2.0/3e8,env.steps+(highest-lowest)/env.dx)
    
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        plt.ylim((env.maxDist+highest-lowest)*2.0/3e8,0)
        draw = np.swapaxes(returnData[j],0,1)
  
        plt.contourf(np.arange(len(files)), ys, draw, 100,norm=MidNorm(np.mean(draw)), cmap="Greys") 
        # normalise data
        #draw = (draw-np.amin(draw))/(np.amax(draw)-np.amin(draw))
        #plt.imshow(draw,aspect="auto")
        #plt.contourf(np.arange(len(files)), ys, draw, 100, cmap="Greys",norm=colors.PowerNorm(1.5))
        #plt.pcolormesh(np.arange(len(files)), ys, draw,cmap="Greys") # alternative cmap is "RdBu_r" where +ve = red
        plt.title(titles[j])
        plt.colorbar()
    if save is not None:
        print "saving"
        plt.savefig(save)
    if display:
        plt.show()
    return returnData

def worker(args):
    global env
    i,name,wave,models,directional,env = args
    with h5py.File(name,"r") as f:
        h = f["meta"][2]
    distance,angle,theta,phi = loadArrays(name)

    ars = np.full((len(models),env.steps),0,float)
    for j in range(len(models)):
        ars[j] = processSlice(distance,angle,theta,phi,models[j],wave,directional=directional)
    return i, h, ars



def wiggle(filename,intensityModel=raySpecular,wave=GaussianDot(),display=True):
    """Calculates the response for a single point file.

    Parameters
    filename - string : Name of the file to generate the response for.
    intensityModel - function (optional) : a function of incidence angle for the
        backscatter from a surface. raySpecular by default i.e. cos(2 theta)
    wave - class instance (optional) : a model of the wave to convolve with the
        response. GaussianDot() by default.
    display - bool (optional) : Whether to plot the data or not. Default is True.

    Returns
    Time series of predicted response if successful, otherwise -1.
    """
    try:
        distance,angle,theta,phi = loadArrays(filename)
    except IOError:
        return fileError(filename)
    ys = processSlice(distance,angle,theta,phi,intensityModel=raySpecular,wave=wave)
    # Clipping
    #ys = np.clip(ys,np.percentile(ys,1),np.percentile(ys,99))
    if display:
        xs = np.linspace(0,env.maxDist*2.0/3e8,env.steps)
        plt.plot(xs,ys,color="black")
        plt.fill_between(xs,ys,0,where=(ys>0),color="black")
        plt.show()
    return ys

def manyWiggle(name,adjusted=False,intensityModel=raySpecular,wave=GaussianDot(),rFactor=0,
               directional = direcNone, compareTo = None):
    '''Plots the response for the points in the given directory side by side.

    Parameters
    name - string : directory to look in for point data. Must contain no other
                    files or directories.
    adjusted - bool : indicates that the responses from each radar point should
                      be aligned vertically. i.e. If the radar took two samples
                      at the same coordinates but different elevations, the
                      response from the surface directly below would be at the
                      same point on the plot for both samples.
    intensityModel - function (optional) : a function of incidence angle for the
        backscatter from a surface. raySpecular by default i.e. cos(2 theta)
    wave - class instance (optional) : a model of the wave to convolve with the
        response. GaussianDot() by default.
    rFactor - float (optional) : How intensity should fall with distance. 0 by default.
    directional - function  (optional) : A function of theta and phi for the
        directivity of the antenna.
    compareTo - 2D float array (optional) : An array of data the same shape as generated
        from the input file. This is plotted behind each estimated response.
    
    Returns
    The 2D array of predicted responses, or the normalised correlation coefficient if
    compareTo was provided. If unsuccessful, returns -1.
    '''
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=figsize)
    try:
        files = listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes "point" prefix
    heights = []

    cells = int(np.ceil(np.sqrt(len(files))))*110
    
    out = np.full((len(files),env.steps),0,float)
    for i in range(len(files)):
        filename = files[i]
        try:
            with h5py.File(name+"/"+filename,"r") as f:
                heights.append(f["meta"][2])
            distance,angle,theta,phi = loadArrays(name+"/"+filename)
        except IOError as e:
            print "Could not read h5py file : "+name+"/"+filename
            return -1
        out[i] = processSlice(distance,angle,theta,phi,intensityModel,wave,rFactor,directional)  

    if adjusted:
        highest = max(heights)
        lowest = min(heights)
        draw = np.full((len(files),env.steps + int((highest-lowest)/env.dx)),0)
        for i in range(len(files)):
            start = int((highest-heights[i])/env.dx)
            draw[i][start:start+env.steps] = out[i]
    else:
        draw = out
        highest, lowest = 0, 0
    ys = np.linspace(0,(env.maxDist+highest-lowest)*2.0/3e8,env.steps+(highest-lowest)/env.dx)
    
    # clipping
    #draw = np.clip(draw,np.percentile(draw,1),np.percentile(draw,99))
    if compareTo is not None:
        compareTo *= np.amax(draw)/np.amax(compareTo)
    
    m = np.amax(abs(draw))*2.0
    for i in range(len(files)):
        plt.plot(draw[i]+m*i,ys,"k-",zorder=5,label="model")
        if compareTo is None:
            plt.fill_betweenx(ys,m*i,m*i+draw[i],where=(draw[i]>0),color='k')
        else:
            plt.plot(compareTo[i]+m*i,ys,"r-",zorder=2,label="reference",alpha=0.3)
            if i == 0:
                plt.legend()
    plt.ylim((env.maxDist+highest-lowest)*2.0/3e8,0)
    plt.xlim(-m,m*len(files))
    plt.gca().axes.get_xaxis().set_visible(False)

    ######## Re-enable plotting once done with tests
    #plt.show()
    ########
    if compareTo is None:
        return draw
    a = draw.reshape(-1)
    b = compareTo.reshape(-1)
    return np.corrcoef(a,b)[0,1]
    
def showWave(wave=GaussianDot()):
    """Plot the given wave over the time the radargram records for."""
    x = np.linspace(-env.maxTime/10.0,env.maxTime,env.steps*2)
    y = wave.amplitude(x)
    plt.plot(x,y)
    plt.plot([wave.align,wave.align],[np.amin(y),np.amax(y)])
    plt.plot([wave.align+env.dt,wave.align+env.dt],[np.amin(y),np.amax(y)])
    plt.show()

def showDirectionality(directional=direcBroad,twoD = False):
    """Creates a spherical plot of the directivity of the antenna.
    Setting twoD causes the model to ignore phi, the azimuth angle."""
    import mpl_toolkits.mplot3d.axes3d as axes3d
    if twoD:
        theta, phi = np.linspace(0, 2 * np.pi, 361), 0
    else:
        theta, phi = np.linspace(0, 2 * np.pi, 61), np.linspace(0, np.pi, 31)
    theta, phi = np.meshgrid(theta,phi)
    r = directional(theta,phi)
    X = r * np.sin(theta) * np.cos(phi)
    Y = r * np.sin(theta) * np.sin(phi)
    Z = r * np.cos(theta)
    fig = plt.figure()
    if twoD:
        ax = fig.add_subplot(1,1,1)
        ax.plot(X[0],Z[0])
        ax.plot([-0.5,0.5],[0,0])
    else:
        ax = fig.add_subplot(1,1,1, projection='3d')
        ax.plot([-0.5,0.5], [0,0], [0,0],linewidth=2.0)
        ax.plot_surface(Z, X, Y, rstride=1, cstride=1,cmap=plt.get_cmap('jet'),antialiased=False,alpha=0.5) 
    plt.show()
    
def fileError(f):
    print "Error in models.py, could not read file/directory: "+f
    print "Please check path/filename entered correctly and data exists"
    return -1
