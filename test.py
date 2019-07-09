import Generate
from os import remove, path
from shutil import rmtree
import numpy as np
from progressbar import progressbar

_startX, _startY = 470305.0, 3094800.0
_endX, _endY = 470038.0, 3095590.0
_steps = 10
_filename = "simple"

# wallPath: (481990, 3085940) to (480892, 3087233)
# lateral: (472971,3083950) to (473882, 3084313)

def setDetailed():
    global _filename, _steps
    _filename = "detailed"
    _steps = 50

def setPath(startX,startY,endX,endY,steps,name):
    global _startX,_startY,_endX,_endY,_steps, _filename
    _startX, _startY = startX, startY
    _endX, _endY = endX, endY
    _steps = steps
    _filename = name

def generateAll():
    xs = np.linspace(_startX,_endX,_steps)
    ys = np.linspace(_startY,_endY,_steps)
    with open("tmp","w") as f:
        f.write(str(_startX)+","+str(_startY)+","+str(_endX)+","+str(_endY)+","+str(_steps)+","+_filename+"\n")
    print "Performing arcpy setup..."
    Generate.Setup()
    print "Complete"
    for i in progressbar(range(_steps)):
        Generate.setPoint(xs[i],ys[i],_filename+"\\point"+str(i))
        Generate.generateMaps()
        with open("tmp","a") as f:
            f.write(str(i)+"\n")
    remove("tmp")
    Generate.finish()

def resume():
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
        if path.exists(_filename+"\\point"+str(done)):
            try:
                rmtree(_filename+"\\point"+str(done))
            except OSError as e:
                print "cleanup failed, please run resume() again"
                return 
        xs = np.linspace(_startX,_endX,_steps)
        ys = np.linspace(_startY,_endY,_steps)
        for i in progressbar(range(done,_steps),min_value=done,max_value=_steps):
            Generate.setPoint(xs[i],ys[i],_filename+"\\point"+str(i))
            Generate.generateMaps()
            with open("tmp","a") as f:
                f.write(str(i)+"\n")
        remove("tmp")
        Generate.finish()
