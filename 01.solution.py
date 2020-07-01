# IMPORTS
import numpy as np
from vpython import *


# CONSTANTS
SIGMA = 0.272
TEMP = 47.0
CELL_SIZE = 1.2 * SIGMA
RANDOM_AREA = 0.8 * SIGMA
CELL_IDENT = (CELL_SIZE - RANDOM_AREA) / 2

GRID_SIZE = 2
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
    """ Строит пустые начальные ускорения """
    return np.empty(N_ATOM*3).reshape(N_ATOM, 3)


def initSpeed():
    """ Строит случайные начальные скорости """

    v_ = np.random.uniform(-1, 1, (N_ATOM, 3))
    # нормировка амплитуд скоростей, чтобы кинетическая энергия соотвествовала начальной температуре
    velsq = np.sum(v_ ** 2)
    aheat = 3.e0*N_ATOM * INTEGRATION_STEP ** 2 * TEMP
    factor = np.sqrt(aheat / velsq)
    
    return v_ * factor


def initForces():
    """ Строит нулевой массив векторов сил """
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
acceleration = initAcceleration()  # массив ускорений
speed = initSpeed()  # массив скоростей
forces = initForces()  # массив векторов сил, действующих на каждую частицу

atoms = createAtomsByPos(coords)
#energy = 0

# LIFECYCLE

for iter in range(ITERATIONS):
    
    scene.waitfor('click')

    
    for i in range(N_ATOM - 1):
        print("coords i",  coords[i])
        xij = coords[i] - coords[i+1:]
        print("xij", xij)
        
        xij[xij < -CELL_SIZE / 2] += CELL_SIZE
        xij[xij > CELL_SIZE / 2] -= CELL_SIZE
        
        rsq = np.sum(xij*xij, axis=1)
        rsqinv = 1.0/rsq
        r6inv = rsqinv*rsqinv*rsqinv

        #enr=4.0*r6inv*(r6inv-1.0)
        force_ = np.einsum('ki,k->ki', xij, rsqinv*48*r6inv*(r6inv-0.5))

        forces[i] += np.sum(force_, axis=0)
        forces[i+1:] -= force_
        #energy=energy+enr
# вычисление (безразмерных) ускорений атомов системы
    acceleration = forces * INTEGRATION_STEP * INTEGRATION_STEP / 2
# вычисление скоростей атомов
    speed += 2.0*acceleration
# новые положения частиц
    coords += speed + acceleration
    coords[coords < 0] += CELL_SIZE
    coords[coords > CELL_SIZE] -= CELL_SIZE

    print(coords[0], speed[0], acceleration[0])
    changePosAtoms(atoms, coords)
# один раз за 5 (kstep) шагов результаты о координатах частиц выводятся в файл,
# и информация о номере шага выводится на экран
# для ускорения расчетов экзаменационного задания вывод в файл можно закомментировать
    
# Обнуляем массив сил
    forces[:, :] = 0
