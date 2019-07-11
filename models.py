import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import listdir
import progressbar
from scipy import signal

_MAXDIST = 30*105.0
_GRANULARITY = 1e-8 # 10ns
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
def Minnart(theta,k): 
    return np.power(np.cos(theta),2*k-2)
def raySpecular(theta,n=1):
    return np.power(np.maximum(0,np.cos(theta*2)),n) 

_wavelength = 50
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
    lim = 1.0/_freq
    def amplitude(self,t):
        return (1-2*np.pi**2*_freq**2*t**2)*np.exp(-np.pi**2*_freq**2*t**2)
class GaussianDot:
    lim = 2.0/(3.0*_freq)
    def amplitude(self,t):
        return -np.exp(0.5)*2*np.pi*_freq*t*np.exp(-2*np.pi**2*_freq**2*t**2)
class Constant:
    lim = _GRANULARITY/2.0
    def amplitude(self,t):
        return 1.0
class IDLWave:
    lim = _MAXTIME
    def __init__(self,c=5.0):
        self.c = float(c)
    def amplitude(self,t):
        return np.exp(-t**2/self.c)

def directionality(theta,phi):
    return np.sin(theta)**2 # IDL model appeared to use sin(theta)**8

def processSlice(filename,intensityModel=raySpecular,wave=GaussianDot()):
    arrays = np.load(filename+"/arrays.npz")
    visible = arrays["visible"]
    distance = arrays["distance"]
    angle = arrays["incidence"]
    theta,phi = None, None
    if "antennaTheta" in list(arrays):
        theta = arrays["antennaTheta"]
        phi = arrays["antennaPhi"]
    height,width = visible.shape[0], visible.shape[1]
    #sample = np.full((_steps),0,"float")
    
    m = (visible == 1) & (distance < _MAXDIST)
    t = (_time(distance[m])/_GRANULARITY).astype(int)
    intensity = intensityModel(angle[m])*directionality(theta[m],phi[m])

    sample = np.array([np.sum(intensity[t==i]) for i in range(_steps)])
    #for i in range(_steps):
    #    sample[i] = np.sum(intensity[t==i])
    low = int(math.floor(-wave.lim/_GRANULARITY))
    high = 1-low
    w = np.full((high-low),0,"float")
    for j in range(high-low):
        w[j] = wave.amplitude((j+low)*_GRANULARITY)
    return np.convolve(sample,w,"same")

#models = [lambertian,lambda x : Minnart(x,2), raySpecular, lambda x : raySpecular(x,2)]
#titles = ["diffuse", "Minnart k=2", "Specular n=1", "Specular n=2"]
models = [raySpecular]
titles = ["Specular n=1"]
#models = [lambda x : Minnart(x,3)]
#titles = ["IDL method"]
def compare(name,adjusted=False,wave=GaussianDot()):
    returnData = []
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(15,10))
    files = listdir(name)
    heights = []
    
    with open(name+"/"+files[0]+"/x_y_z_elevation","r") as f:
        meta = f.read().split(",")
        x0 = float(meta[0])
        y0 = float(meta[1])
    cells = int(np.ceil(np.sqrt(len(models))))*110
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        out = np.full((len(files),_steps),0)
        for i in progressbar.progressbar(range(len(files))):
            filename = files[i]
            with open(name+"/"+filename+"/x_y_z_elevation","r") as f:
                meta = f.read().split(",")
                heights.append(float(meta[2]))
            out[i] = processSlice(name+"/"+filename,models[j],wave) 

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
        plt.ylim((_MAXDIST+highest-lowest)*2.0/3e8,0)

        draw = np.swapaxes(draw,0,1)
        returnData.append(draw)
        plt.contourf(np.arange(len(files)), ys, draw, 100,norm=MidNorm(0), cmap="Greys") 
        # normalise data
        #draw = (draw-np.amin(draw))/(np.amax(draw)-np.amin(draw))
        #plt.imshow(draw,aspect="auto")
        #plt.contourf(np.arange(len(files)), ys, draw, 100, cmap="Greys",norm=colors.PowerNorm(1.5))
        #plt.pcolormesh(np.arange(len(files)), ys, draw,cmap="Greys") # alternative cmap is "RdBu_r" where +ve = red
        plt.title(titles[j])
        plt.colorbar()
    plt.show()
    return returnData

def wiggle(filename,intensityModel=raySpecular,wave=GaussianDot()):
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(15,10))
    ys = processSlice(filename,intensityModel=raySpecular,wave=wave)
    xs = np.linspace(0,_MAXDIST*2.0/3e8,_steps)
    plt.plot(xs,ys,color="black")
    plt.fill_between(xs,ys,0,where=(ys>0),color="black")
    plt.show()

def manyWiggle(name,adjusted=False,intensityModel=raySpecular,wave=GaussianDot()):
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(20,10))
    files = listdir(name)
    heights = []
    
    with open(name+"/"+files[0]+"/x_y_z_elevation","r") as f:
        meta = f.read().split(",")
        x0 = float(meta[0])
        y0 = float(meta[1])
    cells = int(np.ceil(np.sqrt(len(files))))*110
    
    out = np.full((len(files),_steps),0)
    for i in progressbar.progressbar(range(len(files))):
        filename = files[i]
        with open(name+"/"+filename+"/x_y_z_elevation","r") as f:
            meta = f.read().split(",")
            heights.append(float(meta[2]))
        out[i] = processSlice(name+"/"+filename,intensityModel,wave)  

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

    m = np.amax(abs(draw))
    for i in range(len(files)):
        plt.plot(draw[i]+m*i,ys,"k-")
        plt.fill_betweenx(ys,m*i,m*i+draw[i],where=(draw[i]>0),color='k')
    plt.ylim((_MAXDIST+highest-lowest)*2.0/3e8,0)
    plt.xlim(-m,m*len(files))
    plt.gca().axes.get_xaxis().set_visible(False)
    print np.amax(abs(draw))
    plt.show()
    
def showWave(wave=GaussianDot()):
    x = np.linspace(-1.5*wave.lim,1.5*wave.lim,100)
    f = np.vectorize(wave.amplitude)
    plt.plot(x,f(x))
    plt.show()
