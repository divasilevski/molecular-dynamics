import numpy as np
from vpython import *

# CONSTANTS
SIGMA = 0.272

CELL_SIZE_A = 0.9 * SIGMA
CELL_SIZE_B = 1.2 * SIGMA
RANDOM_AREA = 0.8 * SIGMA
CELL_IDENT = (CELL_SIZE_B - RANDOM_AREA) / 2

GRID_SIZE = 10
N_ATOM = GRID_SIZE ** 3

ATOM_RADIUS = CELL_SIZE_B / 10
CAMERA_SIZE = 500

# FUNCTIONS
def calcPosByCellA(cellPos):
    """ Рассчитывает координаты атома по положению ячейки, вариант A """

    initCoord = CELL_SIZE_A * SIGMA * cellPos
    identArray = CELL_SIZE_A * SIGMA / 2 * np.ones(3)

    return initCoord + identArray


def calcPosByCellB(cellPos):
    """ Рассчитывает координаты атома по положению ячейки, вариант B """

    initCoord = CELL_SIZE_B * cellPos
    identArray = CELL_IDENT * np.ones(3)
    randomArray = np.random.uniform(0, RANDOM_AREA, 3)

    return initCoord + identArray + randomArray


def initPositions(calcFunc):
    """ Строит первоначальное расположение атомов """
    initPos = np.empty(N_ATOM * 3).reshape(N_ATOM, 3)  # Выделяем память

    t = 0
    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            for z in range(GRID_SIZE):
                initPos[t] = calcFunc(np.array([x, y, z]))

                t += 1

    return initPos


# VPYTHON FUNCTIONS
def createAtomsByPos(pos):
    """ Создание сфер-атомов """
    points = []
    for i in range(len(pos)):
        P = vector(pos[i][0], pos[i][1], pos[i][2])
        points.append(sphere(pos=P, radius=ATOM_RADIUS))

    return points


def changePosAtoms(atoms, positions):
    """ Изменение позиций сфер-атомов """
    for i in range(len(atoms)):
        P = positions[i]
        atoms[i].pos = vector(P[0], P[1], P[2])


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
scene = canvas(width=CAMERA_SIZE, height=CAMERA_SIZE)
#scene.camera.pos = vector(CAMERA_POS, CAMERA_POS, 0)
scene.autoscale = True

coords = initPositions(calcPosByCellA)
atoms = createAtomsByPos(coords)

scene.waitfor('click')

positions = initPositions(calcPosByCellB)
changePosAtoms(atoms, positions)
