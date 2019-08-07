"""Uses a 'maps.hdf5' file and directory of .hdf5 files for points to model
the radargram/wiggle plots seen. Also contains a range of methods to allow
the model to be altered, such as the backscatter model or the wave to convolve
with the result."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from scipy import signal
import multiprocessing as mp
from progress import progress
import h5py
import align
import modelling

__all__ = ["loadArrays", "processSlice", "models", "titles", "compare", "wiggle", "manyWiggle",
           "showWave", "showDirectionality"]

env = modelling.parameters.env

def _time(distance):
    """Get the time for a wave to return to the radar after reflecting off a
    point 'distance' away."""
    return distance*2.0/3e8


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
        # interpolate so that given value is always 50% and linear either side
        x,y = [self.vmin,self.midpoint,self.vmax],[0,0.5,1]
        return np.ma.masked_array(np.interp(value,x,y),np.isnan(value))

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
def processSlice(distance,angle,theta,phi,intensityModel=modelling.backscatter.raySpecular,
                 wave=modelling.waves.GaussianDot(),rFactor=0.0,directional=modelling.directivity.constant):
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
    t = (_time(distance[m])/env.dt).astype(int) # indices into radargram timesteps
    
    intensity = intensityModel(angle[m]) / np.power(distance[m],rFactor)
    
    # directionality 
    if theta is not None:
        intensity *= directional(theta[m],phi[m])

    # total intensity recieved at each timestep
    sample = np.array([np.sum(intensity[t==i]) for i in range(env.steps)])
    #sample = np.array([np.sum(t==i) for i in range(_steps)])
    #sample = np.array([np.sum(intensity[t==i])/(1+np.sum(t==i)) for i in range(_steps)])
    
    # convolution
    w = wave.amplitude(np.linspace(0.0,env.maxTime,env.steps))
    idx = int(wave.offset/env.dt)
    convol = np.convolve(sample,w)[idx:idx+env.steps]
    
    # Filtering
    #convol = bandpass_filter(convol,0.05,0.4,1.0)

    return convol

# default models to compare radargrams of
#models = [modelling.backscatter.rayModel,modelling.backscatter.ray2Model,
   #       modelling.backscatter.Min2,modelling.backscatter.specular8Model]
#titles = ["Ray tracing n=1","Ray tracing n=2","Minnaert k=2", "spec8"]
models = [modelling.backscatter.Min2]
titles = ["title"]

def compare(name,adjusted=False,wave=modelling.waves.GaussianDot(),save=None,models=models,
            titles=titles,directional=modelling.directivity.constant,display=True, clip = 0
            ,rFactor = 0.0, parallel = True):
    '''Plots the radargram for the points in the given directory once for each
    model given in the models list. Adjusted aligns the y-axis by elevation
    rather than timing. The wave used by the model can also be changed and setting
    "save" to a string means the plot will be saved to that location.

    Parameters
    name - string : directory to look in for point data. Must contain no other
                    files or directories.
    adjusted - bool : The data for each point is shifted to align the response from
                    the surface directly beneath the radar with the top of the plot.
    wave - Wave instance : the wave to convolve over the response.
    save - string (optional) : the name of the file to save the radargram in.
    models - function list (optional) : Functions of incidence angle to use for modelling the response.
    titles - string list (optional) : Names to display above each plot, corresponding to the model
        in the same index in the models list. Both use default lists defined in this module.
    directional - function (optional) : A function for the directionality of the antenna
        in terms of spherical coordinates from the ends of the antenna.
    display - bool (optional) : Display the plot once draw. Default True.
    clip - float (optional) : Clip X percent from either extreme. Default 0.

    Returns
    returnData - 2D float array : The generated radargram data.
    Returns -1 if the method fails.
    '''
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=modelling.parameters.figsize)
    try:
        files = os.listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes 'pointX.hdf5'
    
    returnData = np.full((len(models),len(files),env.steps),0,float) # 3D - many subplots
    
    p = mp.Pool(mp.cpu_count())
    # data needed to run across several processors
    data = [(i,name+"/"+files[i],wave,models,directional,env,rFactor) for i in range(len(files))]
    try:
        if parallel:
            for i, ars in progress(p.imap_unordered(worker,data),len(files)):
                returnData[:,i] = ars
        else: # non-parallel option. Runs on same thread so get proper error reporting
            for j in progress(range(len(files))):
                i,ars = worker(data[j])
                returnData[:,i] = ars
    except IOError as e: 
        p.close()
        print "\nError reading hdf5 file :\n"+e.message
        return -1
    p.close()

    if adjusted: # align first responses with top of radargram
        for i in range(len(models)):
            returnData[i] =  align.minAlign(returnData[i], env.dx, 200.0)

    cells = int(np.ceil(np.sqrt(len(models))))*110
    ys = np.linspace(0, env.maxTime, env.steps)

    # clipping
    if clip:
        returnData = np.clip(returnData,np.percentile(returnData,clip),np.percentile(returnData,100-clip))
    
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        plt.ylim(env.maxTime, 0)
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
    i,name,wave,models,directional,env,rFactor = args
    distance,angle,theta,phi = loadArrays(name)
    ars = np.full((len(models),env.steps),0,float)
    for j in range(len(models)):
        ars[j] = processSlice(distance,angle,theta,phi,models[j],wave,directional=directional,
                              rFactor = rFactor)
    return i, ars

def wiggle(filename,intensityModel=modelling.backscatter.raySpecular,
           wave=modelling.waves.GaussianDot(),display=True):
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
        xs = np.linspace(0,env.maxTime,env.steps)
        plt.plot(xs,ys,color="black")
        plt.fill_between(xs,ys,0,where=(ys>0),color="black")
        plt.show()
    return ys

def manyWiggle(name,adjusted=False,intensityModel=modelling.backscatter.raySpecular,
               wave=modelling.waves.GaussianDot(),rFactor=0,
               directional = modelling.directivity.constant, compareTo = None):
    '''Plots the response for the points in the given directory side by side.

    Parameters
    name - string : directory to look in for point data. Must contain no other
                    files or directories.
    adjusted - bool (optional) : the data is shifted for each point to align the response
                    from the surface directly beneath the radar with the top of the plot.
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
    plt.figure(figsize=modelling.parameters.figsize)
    try:
        files = os.listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes "point" prefix
    heights = []

    cells = int(np.ceil(np.sqrt(len(files))))*110
    
    draw = np.full((len(files),env.steps),0,float)
    for i in range(len(files)):
        filename = files[i]
        try:
            distance,angle,theta,phi = loadArrays(name+"/"+filename)
        except IOError as e:
            print "Could not read h5py file : "+name+"/"+filename
            return -1
        draw[i] = processSlice(distance,angle,theta,phi,intensityModel,wave,rFactor,directional)  

    if adjusted:
        draw =  align.minAlign(draw, env.dx, 200.0)
    ys = np.linspace(0, env.maxTime, env.steps)
    
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
    plt.ylim(env.maxTime,0)
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
    
def showWave(wave=modelling.waves.GaussianDot()):
    """Plot the given wave over the time the radargram records for."""
    x = np.linspace(-env.maxTime/10.0,env.maxTime,env.steps*2)
    y = wave.amplitude(x)
    plt.plot(x,y)
    plt.plot([wave.offset,wave.offset],[np.amin(y),np.amax(y)])
    plt.plot([wave.offset+env.dt,wave.offset+env.dt],[np.amin(y),np.amax(y)])
    plt.show()
    return x,y

def showDirectivity(directional=modelling.directivity.broad, twoD = False):
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
