from vpython import *
scene = canvas(width=500, height=500)
scene.autoscale = True

# CONSTANTS
SIGMA = 0.272
CELL_SIZE = 1.2 * SIGMA
ATOM_RADIUS = CELL_SIZE / 10

data = open("data.coordinates.txt", "r")


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
        points.append(sphere(pos=p, radius=ATOM_RADIUS))
    return points


def changePosition(points, positions):
    for i in range(len(points)):
        points[i].pos = positions[i]


points = createPoints(getNextPosition())
scene.waitfor('click')
while True:
    # sleep(0.5)

    gnp = getNextPosition()
    if gnp:
        changePosition(points, gnp)
    else:
        break

data.close()
