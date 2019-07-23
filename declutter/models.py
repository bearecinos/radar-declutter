import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import listdir
from scipy import signal
import multiprocessing as mp
from progress import progress

# plt.ion makes plots interactive (nonblocking) but may harm performance
# and/or stability

_MAXDIST = 3000.0 
_GRANULARITY = 1e-8 # 10ns
_SPACE_GRANULARITY = _GRANULARITY*3e8
_MAXTIME = 2*_MAXDIST/3e8

_steps = int(_MAXTIME/_GRANULARITY)

_MetaFilePath = "/x_y_z"

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
    # exp( - (findgen(echo_size)-float(i))^2/50. )
    # timestep was 240us (2.4e-4 seconds echo length) / 2048 granularity
    # = 1.172e-7 = 117ns - much longer than our sampling rate
    # scaling = 1.172e-7
    lim = _MAXTIME
    def __init__(self,c=50.0):
        self.c = float(c)
    def amplitude(self,t):
        return np.exp(-(t)**2/self.c)
class RC:
    def __init__(self,c=5.0):
        self.c = float(c)
    lim = 30e-6
    def amplitude(self,t):
        if t < 0:
            return 0.0
        return np.exp(-t/(self.c*1e-6))

# Squared to account for also being receiver?
def directionality(theta,phi):
    # IDL model appeared to use sin(theta)**8
    # sin(3*theta) term adds side lobes of lesser intensity    
    #return np.sin(theta)**2*np.sin(3*theta)**2
    #return np.sin(theta)**4
    return abs(np.sin(theta))**0.5

# Replace GaussianDot with RC
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
    convol = np.convolve(sample,w) # no longer 'same' mode, now full so must crop
    if len(convol) > _steps:
        return convol[len(w)/2:len(w)/2+_steps]
    #return np.convolve(sample,w,"same") # TODO: enforce cropping to length of sample
    return convol

#models = [lambertian,lambda x : Minnaert(x,2), raySpecular, lambda x : raySpecular(x,2)]
#titles = ["diffuse", "Minnaert k=2", "Specular n=1", "Specular n=2"]
models = [raySpecular]
titles = ["Specular n=1"]
#models = [lambda x : Minnart(x,3)] 
#titles = ["IDL method"]
def compare(name,adjusted=False,wave=GaussianDot()):
    returnData = []
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(12,8))
    try:
        files = listdir(name)
    except OSError as e:
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:])) # assumes 'point' prefix
    heights = []
    
    with open(name+"/"+files[0]+_MetaFilePath,"r") as f:
        meta = f.read().split(",")
        x0 = float(meta[0])
        y0 = float(meta[1])
    cells = int(np.ceil(np.sqrt(len(models))))*110
    for j in range(len(models)):
        plt.subplot(cells+j+1)
        out = np.full((len(files),_steps),0,float)

        # Bit to parallelize
        # Note - must pass actual function, not a lambda expression
        p = mp.Pool(mp.cpu_count())
        data = [(i,(name+"/"+files[i],models[j],wave)) for i in range(len(files))] 
        for i,h,ar in progress(p.imap_unordered(worker,data),len(files)):
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
        # TODO: no longer want MidNorm? # norm=MidNorm(0)
        # mean is about 10^-10 for GaussianDot anyway
        plt.contourf(np.arange(len(files)), ys, draw, 100,norm=MidNorm(np.mean(draw)), cmap="Greys") 
        # normalise data
        #draw = (draw-np.amin(draw))/(np.amax(draw)-np.amin(draw))
        #plt.imshow(draw,aspect="auto")
        #plt.contourf(np.arange(len(files)), ys, draw, 100, cmap="Greys",norm=colors.PowerNorm(1.5))
        #plt.pcolormesh(np.arange(len(files)), ys, draw,cmap="Greys") # alternative cmap is "RdBu_r" where +ve = red
        plt.title(titles[j])
        plt.colorbar()
        print "Plot complete: "+titles[j]
    plt.show()
    return returnData

def worker(args):
    i = args[0]
    args = args[1]
    with open(args[0]+_MetaFilePath,"r") as f:
        h = float(f.read().split(",")[2])
    return i, h, processSlice(*args) 

def wiggle(filename,intensityModel=raySpecular,wave=GaussianDot()):
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=(12,8))
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
        return fileError(e.filename)
    files.sort(key=lambda x : int(x[5:])) # assumes "point" prefix
    heights = []
    
    with open(name+"/"+files[0]+"/x_y_z","r") as f:
        meta = f.read().split(",")
        x0 = float(meta[0])
        y0 = float(meta[1])
    cells = int(np.ceil(np.sqrt(len(files))))*110
    
    out = np.full((len(files),_steps),0)
    for i in range(len(files)):
        filename = files[i]
        #filename = "point"+str(i)
        with open(name+"/"+filename+_MetaFilePath,"r") as f:
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
    plt.show()
    
def showWave(wave=GaussianDot()):
    x = np.linspace(-1.5*wave.lim,1.5*wave.lim,500)
    f = np.vectorize(wave.amplitude)
    y = f(x)
    plt.plot(x,y)
    plt.plot([0,0],[np.amin(y),np.amax(y)])
    plt.plot([_GRANULARITY,_GRANULARITY],[np.amin(y),np.amax(y)])
    plt.show()
def showDirectionality():
    import mpl_toolkits.mplot3d.axes3d as axes3d
    theta, phi = np.linspace(0, 2 * np.pi, 30), np.linspace(0, np.pi, 15)
    theta, phi = np.meshgrid(theta,phi)
    r = directionality(theta,phi)
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
