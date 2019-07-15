import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import listdir
from scipy import signal
import multiprocessing as mp

# plt.ion makes plots interactive (nonblocking) but may harm performance
# and/or stability

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
def Minnaert(theta,k): 
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

# Squared to account for also being receiver?
def directionality(theta,phi):
    # IDL model appeared to use sin(theta)**8
    # sin(3*theta) term adds side lobes of lesser intensity    
    #return np.sin(theta)**2*np.sin(3*theta)**2
    return np.sin(phi)**4

def processSlice(filename,intensityModel=raySpecular,wave=GaussianDot()):
    try:
        arrays = np.load(filename+"/arrays.npz")
    except IOError:
        print "Error in models.py, could not load file "+filename+"/arrays.npz"
        print "Note, if trying to process a directory of points, must be no other files/directories in directory."
    visible = arrays["visible"]
    distance = arrays["distance"]
    angle = arrays["incidence"]
    theta,phi = None, None
    if "antennaTheta" in list(arrays):
        theta = arrays["antennaTheta"]
        phi = arrays["antennaPhi"]
    height,width = visible.shape[0], visible.shape[1]
        
    m = (visible == 1) & (distance < _MAXDIST)
    t = (_time(distance[m])/_GRANULARITY).astype(int)
    intensity = intensityModel(angle[m])
    if theta is not None:
        intensity *= directionality(theta[m],phi[m])

    sample = np.array([np.sum(intensity[t==i]) for i in range(_steps)])
    
    low = int(math.floor(-wave.lim/_GRANULARITY))
    high = 1-low
    w = np.full((high-low),0,"float")
    for j in range(high-low):
        w[j] = wave.amplitude((j+low)*_GRANULARITY)
    return np.convolve(sample,w,"same")

#models = [lambertian,lambda x : Minnaert(x,2), raySpecular, lambda x : raySpecular(x,2)]
#titles = ["diffuse", "Minnaert k=2", "Specular n=1", "Specular n=2"]
models = [raySpecular]
titles = ["Specular n=1"]
#models = [lambda x : Minnart(x,3)] 
#titles = ["IDL method"]
def compare(name,adjusted=False,wave=GaussianDot()):
    returnData = []
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(15,10))
    try:
        files = listdir(name)
    except OSError as e:
        fileError(e.filename)
    heights = []
    
    with open(name+"/"+files[0]+"/x_y_z_elevation","r") as f:
        meta = f.read().split(",")
        x0 = float(meta[0])
        y0 = float(meta[1])
        cells = int(np.ceil(np.sqrt(len(models))))*110
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        out = np.full((len(files),_steps),0)

        # Bit to parallelize
        p = mp.Pool(mp.cpu_count())
        data = [(i,(name+"/"+files[i],models[j],wave)) for i in range(len(files))]
        for i,h,ar in p.imap_unordered(worker,data):
            heights.append(h)
            out[i] = ar

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

def worker(args):
    i = args[0]
    args = args[1]
    with open(args[0]+"/x_y_z_elevation","r") as f:
        h = float(f.read().split(",")[2])
    return i, h, processSlice(*args) 

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
    try:
        files = listdir(name)
    except OSError as e:
        fileError(e.filename)
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
def showDirectionality():
    import mpl_toolkits.mplot3d.axes3d as axes3d
    theta, phi = np.linspace(0, 2 * np.pi, 90), np.linspace(0, np.pi, 45)
    theta, phi = np.meshgrid(theta,phi)
    r = directionality(theta,phi)
    X = r * np.sin(phi) * np.cos(theta)
    Y = r * np.sin(phi) * np.sin(theta)
    Z = r * np.cos(phi)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot([0,0], [0,0], [-0.4,0.4],linewidth=2.0)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,antialiased=False,cmap=plt.get_cmap('jet'),alpha=0.5)
    plt.show()

def fileError(f):
    print "Error in models.py, could not find file: "+f
    print "Please check path/filename entered correct and data exists"
    return -1
