# IMPORTS
import numpy as np
from vpython import *


# CONSTANTS Вариант 5 (Ne)
GRID_SIZE = 10
SIGMA = 0.272
N_ATOM = GRID_SIZE ** 3

ATOM_RADIUS = 1.2 * SIGMA / 10


# FUNCTIONS
def calcPosByCell(cellPos):
    """ Рассчитывает координаты атома по положению ячейки """

    initCoord = 1.2 * SIGMA * cellPos
    identArray = 0.2 * SIGMA * np.ones(3)
    randomArray = np.random.uniform(0, 0.8 * SIGMA, 3)

    return initCoord + identArray + randomArray


def initPositions():
    """ Строит первоначальное расположение атомов """
    initPos = np.empty(N_ATOM * 3).reshape(N_ATOM, 3)  # Выделяем память

    t = 0
    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            for z in range(GRID_SIZE):
                initPos[t] = calcPosByCell(np.array([x, y, z]))
                t += 1

    return initPos


# VPYTHON FUNCTIONS
def createPoints(positions):
    points = []
    for p in positions:
        p = vector(p[0], p[1], p[2])
        points.append(sphere(pos=p, radius=ATOM_RADIUS))
    return points


def changePosition(points, positions):
    for i in range(len(points)):
        p = positions[i]
        p = vector(p[0], p[1], p[2])
        points[i].pos = p


# LIFECYCLE
scene = canvas(width=500, height=500)
scene.camera.pos = vector(3, 3, 3)
scene.autoscale = False

points = createPoints(initPositions())
