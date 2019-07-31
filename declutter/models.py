import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import listdir
from scipy import signal
import multiprocessing as mp
from progress import progress
import h5py

# plt.ion makes plots interactive (nonblocking) but may harm performance
# and/or stability

# Dist assuming 3x10^8. If basing on radargram, account for 1.67x10^8
# i.e. 1km at ice speed is actually about 1.8km away

#_MAXDIST = 3000.0
_MAXDIST = 701 * 1.25e-8 * 1.5e8 # from sample radargram
_GRANULARITY = 1.25e-8 # 10ns
_SPACE_GRANULARITY = _GRANULARITY*3e8
_MAXTIME = 2*_MAXDIST/3e8

_steps = int(_MAXTIME/_GRANULARITY)

def _time(distance):
    """Get the time for a wave to return to the radar after reflecting off a
    point 'distance' away."""
    return distance*2.0/3e8

# surface models
def lambertian(angle):
    return np.cos(angle)
def lambSpec(angle): # specular but adding fact that area seen is like cos theta
    return lambertian(angle)*raySpecular(angle)
def Henyey_Green(g):
    # backscatter so theta = pi, cos(pi) = -1
    # 1 + g^2 + 2*g = (1+g)^2
    # just constant per surface as incidence angle irrelevant
    return (1-g*g)*np.power(1+g,-3.0)
def Minnaert(theta,k): 
    return np.power(np.cos(theta),2*k-2)
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

########## UNDO?
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

# wave models
class Ricker:
    delt = 1.0/_freq # time of centre point of wave
    align = 0.0 # time that should align to pulse at t=0
    def amplitude(self,t):
        return ((1-2*np.pi**2*_freq**2*(t-self.delt)**2)*
                np.exp(-np.pi**2*_freq**2*(t-self.delt)**2))
class GaussianDot:
    delt = 2.0/(3.0*_freq)
    align = 0.0
    def amplitude(self,t):
        return (-np.exp(0.5)*2*np.pi*_freq*(t-self.delt)*
                np.exp(-2*np.pi**2*_freq**2*(t-self.delt)**2))
class Constant:
    delt = 0.0
    align = 0.0
    def amplitude(self,t):
        return t < _GRANULARITY 
class IDLWave:
    # exp( - (findgen(echo_size)-float(i))^2/50. )
    # timestep was 240us (2.4e-4 seconds echo length) / 2048 granularity
    # = 1.172e-7 = 117ns - much longer than our sampling rate
    # scaling = 1.172e-7
    align = _MAXTIME / 2.0
    delt = align
    def __init__(self,c=50.0):
        self.c = float(c)
    def amplitude(self,t):
        return np.exp(-(t-self.delt)**2/self.c)
class RC:
    align = 0.0
    def __init__(self,c=5.0):
        self.c = float(c)
    def amplitude(self,t):
        if t < 0:
            return 0.0
        return np.exp(-t/(self.c*1e-6))


class Other:
    delt = 1.0/_freq
    align = 0.0
    def amplitude(self,t):
        p = 2*np.pi*_freq*(t-self.delt)
        return (p**4-10*p**2+8)*np.exp(-p*p/4.0)/8.0
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

# Squared to account for also being receiver?
def directionality(theta,phi):
    # IDL model appeared to use sin(theta)**8
    # sin(3*theta) term adds side lobes of lesser intensity    
    #return np.sin(theta)**2*np.sin(3*theta)**2
    #return np.sin(theta)**4
    return abs(np.sin(theta))**0.5

def direcNone(theta,phi):
    return 1
def direcBroad(theta,phi):
    return abs(np.sin(theta))**0.5
def direcIDL(theta,phi):
    return np.sin(theta)**2
def direcLobes(theta,phi):
    return np.sin(theta)**2*np.sin(3*theta)**2

def loadArrays(filename):
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
    
    m = (distance < _MAXDIST)
    t = (_time(distance[m])/_GRANULARITY).astype(int)
    # Modify to reduce as distance increases
    # intensity = intensityModel(angle[m]) 
    intensity = intensityModel(angle[m]) / np.power(distance[m],rFactor)
    
    # directionality or not? ##################
    # input function to use
    if theta is not None:
        intensity *= directional(theta[m],phi[m])
    
    ###########################################

    # TODO: Decide if this is best way to do wave
    # TODO: Add in filtering or not?

    ##### UNDO - normalizing by number of points in bin #######
    sample = np.array([np.sum(intensity[t==i]) for i in range(_steps)])
    #sample = np.array([np.sum(t==i) for i in range(_steps)])
    #sample = np.array([np.sum(intensity[t==i])/(1+np.sum(t==i)) for i in range(_steps)])
    ###########################################################


    # CONVOLUTION
    w = wave.amplitude(np.linspace(0.0,_MAXTIME,_steps))
    idx = int(wave.align/_GRANULARITY)
    convol = np.convolve(sample,w)[idx:idx+_steps]
    
    # Filtering
    #convol = bandpass_filter(convol,0.05,0.4,1.0)

    return convol

#models = [lambertian,lambda x : Minnaert(x,2), raySpecular, lambda x : raySpecular(x,2)]
#titles = ["diffuse", "Minnaert k=2", "Specular n=1", "Specular n=2"]
models = [rayModel,ray2Model,specular2Model,specular4Model]
titles = ["Ray tracing n=1","Ray tracing n=2","sin(theta)^2","sin(theta)^4"]

def compare(name,adjusted=False,wave=GaussianDot(),save=None):
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
    '''
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(12,8))
    try:
        files = listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes 'pointX.hdf5'
    heights = []
    returnData = np.full((len(models),len(files),_steps),0,float) # 3D - many plots
    
    p = mp.Pool(mp.cpu_count())
    data = [(i,name+"/"+files[i],wave,models) for i in range(len(files))]

    for i, h, ars in progress(p.imap_unordered(worker,data),len(files)):
        # need to handle case where result is invalid: i = -1
        returnData[:,i] = ars
        heights.append(h)
    p.close()

    highest,lowest = 0,0
    if adjusted:
        highest = max(heights)
        lowest = min(heights)
        draw = np.full((len(models),n,_steps + int((highest-lowest)/_SPACE_GRANULARITY)),0)
        for i in range(len(files)):
            start = int((highest-heights[i])/_SPACE_GRANULARITY)
            draw[:,i,start:start+_steps] = returnData[:,i]
        returnData = draw

    cells = int(np.ceil(np.sqrt(len(models))))*110
    ys = np.linspace(0,(_MAXDIST+highest-lowest)*2.0/3e8,_steps+(highest-lowest)/_SPACE_GRANULARITY)
    
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        plt.ylim((_MAXDIST+highest-lowest)*2.0/3e8,0)
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
    plt.show()
    return returnData

def worker(args):
    i,name,wave,models = args
    with h5py.File(name,"r") as f:
        h = f["meta"][2]
    distance,angle,theta,phi = loadArrays(name)

    ars = np.full((len(models),_steps),0,float)
    for j in range(len(models)):
        ars[j] = processSlice(distance,angle,theta,phi,models[j])
    return i, h, ars



def wiggle(filename,intensityModel=raySpecular,wave=GaussianDot()):
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(12,8))
    distance,angle,theta,phi = loadArrays(filename)
    ys = processSlice(distance,angle,theta,phi,intensityModel=raySpecular,wave=wave)
    # Clipping
    #ys = np.clip(ys,np.percentile(ys,1),np.percentile(ys,99))
    xs = np.linspace(0,_MAXDIST*2.0/3e8,_steps)
    plt.plot(xs,ys,color="black")
    plt.fill_between(xs,ys,0,where=(ys>0),color="black")
    plt.show()

def manyWiggle(name,adjusted=False,intensityModel=raySpecular,wave=GaussianDot(),rFactor=0,
               directional = direcNone, compareTo = None):
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(10.5,5.5))
    try:
        files = listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:-5])) # assumes "point" prefix
    heights = []

    cells = int(np.ceil(np.sqrt(len(files))))*110
    
    out = np.full((len(files),_steps),0,float)
    for i in range(len(files)):
        filename = files[i]
        with h5py.File(name+"/"+filename,"r") as f:
            heights.append(f["meta"][2])
        distance,angle,theta,phi = loadArrays(name+"/"+filename)
        out[i] = processSlice(distance,angle,theta,phi,intensityModel,wave,rFactor,directional)  

    if adjusted:
        highest = max(heights)
        lowest = min(heights)
        draw = np.full((len(files),_steps + int((highest-lowest)/_SPACE_GRANULARITY)),0)
        for i in range(len(files)):
            start = int((highest-heights[i])/_SPACE_GRANULARITY)
            draw[i][start:start+_steps] = out[i]
    else:
        draw = out
        highest, lowest = 0, 0
    ys = np.linspace(0,(_MAXDIST+highest-lowest)*2.0/3e8,_steps+(highest-lowest)/_SPACE_GRANULARITY)
    
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
    plt.ylim((_MAXDIST+highest-lowest)*2.0/3e8,0)
    plt.xlim(-m,m*len(files))
    plt.gca().axes.get_xaxis().set_visible(False)

    ######## Re-enable plotting once done with tests
    #plt.show()
    ########
    a = draw.reshape(-1)
    b = compareTo.reshape(-1)
    return np.corrcoef(a,b)[0,1]
    
def showWave(wave=GaussianDot()):
    x = np.linspace(-_MAXTIME/10.0,_MAXTIME,_steps)
    y = wave.amplitude(x)
    
    #f = np.vectorize(wave.amplitude)
    #y = f(x)
    plt.plot(x,y)
    plt.plot([wave.align,wave.align],[np.amin(y),np.amax(y)])
    plt.plot([wave.align+_GRANULARITY,wave.align+_GRANULARITY],[np.amin(y),np.amax(y)])
    plt.show()
def showDirectionality(directional=direcBroad):
    import mpl_toolkits.mplot3d.axes3d as axes3d
    theta, phi = np.linspace(0, 2 * np.pi, 30), np.linspace(0, np.pi, 15)
    theta, phi = np.meshgrid(theta,phi)
    r = directional(theta,phi)
    X = r * np.sin(theta) * np.cos(phi)
    Y = r * np.sin(theta) * np.sin(phi)
    Z = r * np.cos(theta)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot([-0.5,0.5], [0,0], [0,0],linewidth=2.0)
    ax.plot_surface(Z, X, Y, rstride=1, cstride=1,cmap=plt.get_cmap('jet'),antialiased=False,alpha=0.5) 
    plt.show()

def fileError(f):
    print "Error in models.py, could not find file: "+f
    print "Please check path/filename entered correct and data exists"
    return -1
