from vpython import *

# CONSTANTS
SIGMA = 0.272
CELL_SIZE = 1.2 * SIGMA
ATOM_RADIUS = CELL_SIZE / 10
FILE = "01"

dataCoordinates = open("data//data.coordinates." + FILE + ".txt", "r")
dataPlot = open("data//data.plot." + FILE + ".txt", "r")


def splitToVector(splitLine):
    x = float(splitLine[0])
    y = float(splitLine[1])
    z = float(splitLine[2])
    return vector(x, y, z)


def getNextPosition():
    positions = []
    line = dataCoordinates.readline()
    while line:
        splitLine = line.split()

        if splitLine:
            positions.append(splitToVector(splitLine))
        else:
            return positions

        line = dataCoordinates.readline()


def createPoints(positions):
    points = []
    for p in positions:
        points.append(sphere(pos=p, radius=ATOM_RADIUS))
    return points


def changePosition(points, positions):
    for i in range(len(points)):
        points[i].pos = positions[i]

scene = canvas(width=500, height=500)
scene.autoscale = True
graph(width = 600, height = 200, title = dataPlot.readline().rstrip())
curve = gcurve(color=color.orange)

points = createPoints(getNextPosition())
scene.waitfor('click')

iter = 0
while True:
    # sleep(0.5)

    gnp = getNextPosition()
    if gnp:
        changePosition(points, gnp)
    else:
        break
    
    curve.plot(pos=(iter, float(dataPlot.readline())))
    iter += 1

dataCoordinates.close()
dataPlot.close()