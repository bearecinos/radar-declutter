from os import remove, path
from shutil import rmtree
import numpy as np
import stateless
import multiprocessing as mp

_startX, _startY = 470305.0, 3094800.0
_endX, _endY = 470038.0, 3095590.0
_steps = 10
_filename = "simple"

# wallPath: (481990, 3085940) to (480892, 3087233)
# lateral: (472971,3083950) to (473882, 3084313)

def setDetailed():
    """Sets the pathname to 'detailed' and number of samples to 50."""
    global _filename, _steps
    _filename = "detailed"
    _steps = 50

def setPath(startX,startY,endX,endY,steps,name):
    """Select the endpoints of a straight path, the resolution of the path, and its name.
    No option to set elevation, 100m above ground at every point.
    Parameters:
    startX float : initial x-coordinate.
    startY float : initial y-coordinate.
    endX float : final x-coordinate.
    endY float : final y-coordiante.
    steps int : number of points to sample.
    name string : name of folder to store data in."""    
    global _startX,_startY,_endX,_endY,_steps, _filename
    _startX, _startY = startX, startY
    _endX, _endY = endX, endY
    _steps = steps
    _filename = name

def generateAll():
    """Generates data for all equally spaced points along the specified straight path."""
    xs = np.linspace(_startX,_endX,_steps)
    ys = np.linspace(_startY,_endY,_steps)
    with open("tmp","w") as f:
        f.write(str(_startX)+","+str(_startY)+","+str(_endX)+","+str(_endY)+","+str(_steps)+","+_filename+"\n")
    stateless.Setup()
    for i in range(_steps):
        stateless.generateMaps(xs[i],ys[i],_filename+"/point"+str(i))
        with open("tmp","a") as f:
            f.write(str(i)+"\n")
    remove("tmp")

def workerCall(args):
    x,y,name = args
    stateless.generateMaps(x,y,name)

def parallelAll():
    xs = np.linspace(_startX,_endX,_steps)
    ys = np.linspace(_startY,_endY,_steps)
    
    pool = mp.Pool(mp.cpu_count())
    data = [(x,y,_filename+"/point"+str(i)) for x,y,i in zip(xs,ys,np.arange(_steps))]
    fail = False
    for r in pool.imap_unordered(workerCall,data):
        if r == -1:
            fail = True
    pool.close()
    if fail:
        print "Failed, setup didn't run correctly. Likely issue with state copying between processes"
        return -1
    return 0   
        

def resume():
    """Read the tmp file created when points are generated and use this to resume a partially complete path.
    Note : This does not check the tmp file first exists so will fail if the most recent path completed successfully."""
    with open("tmp","r") as f:
        lines = f.read().split("\n")
    done = len(lines)-2
    meta = lines[0].split(",")
    _startX, _startY = float(meta[0]), float(meta[1])
    _endX, _endY = float(meta[2]), float(meta[3])
    _steps, _filename = int(meta[4]), meta[5]
    if done <= 0:    
        generateAll()
    else:
        if path.exists(_filename+"/point"+str(done)):
            try:
                rmtree(_filename+"/point"+str(done))
            except OSError as e:
                print "cleanup failed, please run resume() again"
                return 
        xs = np.linspace(_startX,_endX,_steps)
        ys = np.linspace(_startY,_endY,_steps)
        stateless.Setup()
        for i in range(done,_steps):
            stateless.generateMaps(xs[i],ys[i],_filename+"/point"+str(i))
            with open("tmp","a") as f:
                f.write(str(i)+"\n")
        remove("tmp")

if __name__=="__main__":
    from time import clock
    setDetailed()
    t0 = clock()
    generateAll()
    result = clock() - t0
    _filename = "detailed1"
    t0 = clock()
    parallelAll()
    result2 = clock() - t0
    print "Serial: "+str(result)
    print "Parallel: "+str(result2)
