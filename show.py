from vpython import *
scene = canvas(width=500, height=500)
scene.camera.pos = vector(3, 3, 3)
scene.autoscale = False

data = open("data.txt", "r")

def splitToVector(splitLine):
    x = float(splitLine[0])
    y = float(splitLine[1])
    z = float(splitLine[2])
    return vector(x, y, z)

def getNextPosition():
    positions = []
    line = data.readline()
    while line:
        splitLine = line.split()
        
        if splitLine: 
            positions.append(splitToVector(splitLine))
        else: 
            return positions
            
        line = data.readline()

def createPoints(positions):
    points = []
    for p in positions:
        points.append(sphere(pos=p, radius=0.1))
    return points

def changePosition(points, positions):
    for i in range(len(points)):
        points[i].pos = positions[i]


points = createPoints(getNextPosition())
scene.waitfor('click')
while True:
    # sleep(0.5)
    
    gnp = getNextPosition()
    if gnp: changePosition(points, gnp)
    else: break

data.close()