# IMPORTS
import numpy as np
from vpython import *


# CONSTANTS
SIGMA = 0.272
CELL_SIZE = 1.2 * SIGMA
RANDOM_AREA = 0.8 * SIGMA
CELL_IDENT = (CELL_SIZE - RANDOM_AREA) / 2

GRID_SIZE = 10
N_ATOM = GRID_SIZE ** 3

ATOM_RADIUS = CELL_SIZE / 10
CAMERA_SIZE = 500
CAMERA_POS = CELL_SIZE * GRID_SIZE / 2

PENTA_TIME = 3
INTEGRATION_STEP = 0.001
ITERATIONS = int(PENTA_TIME / INTEGRATION_STEP)


# FUNCTIONS
def calcPosByCell(cellPos):
    """ Рассчитывает координаты атома по положению ячейки """

    initCoord = CELL_SIZE * cellPos
    identArray = CELL_IDENT * np.ones(3)
    randomArray = np.random.uniform(0, RANDOM_AREA, 3)

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


def initAcceleration():
    """ Строит начальные ускорения """
    return np.empty(N_ATOM*3).reshape(N_ATOM, 3)


def initForces():
    """ Строит начальный массив векторов сил """
    return np.zeros(N_ATOM*3).reshape(N_ATOM, 3)


# VPYTHON FUNCTIONS
def createAtomsByPos(positions):
    """ Создание сфер-атомов """
    points = []
    for P in positions:
        P = vector(P[0], P[1], P[2])
        points.append(sphere(pos=P, radius=ATOM_RADIUS))
    return points


def changePosAtoms(atoms, positions):
    """ Изменение позиций сфер-атомов """
    for i in range(len(atoms)):
        P = positions[i]
        atoms[i].pos = vector(P[0], P[1], P[2])


# PREPARATION
scene = canvas(width=CAMERA_SIZE, height=CAMERA_SIZE)
scene.camera.pos = vector(CAMERA_POS, CAMERA_POS, 0)

coords = initPositions()  # массив радиус-векторов частиц
boost = initAcceleration()  # массив ускорений
forces = initForces()  # массив векторов сил, действующих на каждую частицу

atoms = createAtomsByPos(coords)

# LIFECYCLE
# scene.waitfor('click')
for i in range(ITERATIONS):
    print(i)
