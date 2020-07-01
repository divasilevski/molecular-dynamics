# IMPORTS
import numpy as np
from vpython import *


# CONSTANTS Вариант 5 (Ne)
GRID_SIZE = 10
SIGMA = 0.272
N_ATOM = GRID_SIZE ** 3


# FUNCTIONS
def calcPosByCell(cellPos):
    """ Рассчитывает координаты атома по положению ячейки"""

    return 1.2*SIGMA*cellPos + 0.6*SIGMA*np.array([1,1,1])


def initPositions():
    """ Строит первоначальное расположение атомов """

    initPos = np.empty(N_ATOM*3).reshape(N_ATOM, 3)  # Выделяем память

    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            for z in range(GRID_SIZE):
                initPos[x+y+z] = calcPosByCell(np.array([x, y, z]))

    return initPos


# VPYTHON
print(initPositions())

# LIFECYCLE
