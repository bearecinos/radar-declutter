import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt

def gaussian(order=2,sigma=1.0):
    n = 2*order+1
    k = np.full((n,n),0)
    for i,j in np.ndindex(n,n):
        k[i][j] = np.exp(-((i-order-1.0)**2+(j-order-1.0)**2)/(2.0*sigma**2))/(2*np.pi*sigma**2)
    return k

def bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return signal.lfilter(b, a, data)

def compass(theta): # 0 = N, 1 = NE, 2 = E, 3 = SE
    return int((180.0*theta/np.pi+22.5)/45) % 4

def findEdges(data,order=2):
    k = gaussian(order)
    filtered = ndimage.convolve(data,k)[:-2*order,:-2*order] # edges without enough data cut off

    Lx = np.array([[1,0,-1],[2,0,-2],[1,0,-1]]) # Sobel operator
    Ly = Lx.swapaxes(0,1)

    Gx = ndimage.convolve(filtered,Lx)[:-2,:-2]
    Gy = ndimage.convolve(filtered,Ly)[:-2,:-2]

    G = np.hypot(Gx,Gy)
    Theta = np.arctan2(Gy,Gx)
    dirs = np.vectorize(compass)(Theta)

    boundary = np.percentile(G,90)
    plt.imshow(G>boundary,aspect="auto")
    plt.show()
    onOff = (G>boundary).astype(float)

def recurse(onOff,groups,v,ls):
    ymax,xmax = onOff.shape()
    ymax -= 1
    xmax -= 1
    while ls != []:
        x,y = ls.pop()
        if x > 0:
            if onOff[y][x-1] and groups[y][x-1] == -1:
                groups[y][x-1] = v
                ls.append((x-1,y))
        if x < xmax:
            if onOff[y][x+1] and groups[y][x+1] == -1:
                groups[y][x+1] = v
                ls.append((x+1,y))
        if y > 0:
            if onOff[y-1][x] and groups[y-1][x] == -1:
                groups[y-1][x] = v
                ls.append((x,y-1))
        if y < ymax:
            if onOff[y+1][x] and groups[y+1][x] == -1:
                groups[y+1][x] = v
                ls.append((x,y+1))

def grouping(onOff):
    ymax,xmax = onOff.shape()
    groups = np.full_like(onOff,-1)
    g = 0
    for x,y in np.ndindex(xmax,ymax):
        if onOff[y][x] and groups[y][x] == -1:
            groups[y][x] = g
            recurse(g,[(x,y)])
            g+=1

order = 2
k = gaussian(order)
data = np.load("radargram.npy")

filtered = ndimage.convolve(data,k)[:-2*order,:-2*order] # edges without enough data cut off

Lx = np.array([[1,0,-1],[2,0,-2],[1,0,-1]]) # Sobel operator
Ly = Lx.swapaxes(0,1)

Gx = ndimage.convolve(filtered,Lx)[:-2,:-2]
Gy = ndimage.convolve(filtered,Ly)[:-2,:-2]

G = np.hypot(Gx,Gy)
Theta = np.arctan2(Gy,Gx)
dirs = np.vectorize(compass)(Theta)
# try low pass filtering whole image?
