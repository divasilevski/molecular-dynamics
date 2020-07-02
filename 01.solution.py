# IMPORTS
import numpy as np
from vpython import *


# CONSTANTS
SIGMA = 0.272
BOILING_POINT = 27.1 # https://en.wikipedia.org/wiki/Neon
MELTING_POINT = 24.55
TEMP_EPS = 47.0
TEMP = BOILING_POINT / TEMP_EPS # температура релаксации
TEMP_0 = MELTING_POINT / TEMP_EPS
TEMP_1 = 2 * TEMP

CELL_SIZE = 1.2 * SIGMA
RANDOM_AREA = 0.8 * SIGMA
CELL_IDENT = (CELL_SIZE - RANDOM_AREA) / 2

GRID_SIZE = 10
COUB_SIZE = CELL_SIZE * GRID_SIZE
N_ATOM = GRID_SIZE ** 3

ATOM_RADIUS = CELL_SIZE / 10
CAMERA_SIZE = 500
CAMERA_POS = CELL_SIZE * GRID_SIZE / 2

FULL_TIME = 3.5e-5
INTEGRATION_STEP = 3.5e-8
ITERATIONS = int(FULL_TIME / INTEGRATION_STEP)


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
    aheat = 3.0 * N_ATOM * INTEGRATION_STEP ** 2 * TEMP
    factor = np.sqrt(aheat / velsq)

    return v_ * factor


def initForces():
    """ Строит нулевой массив векторов сил """
    return np.zeros(N_ATOM*3).reshape(N_ATOM, 3)


def calculation(coords, energy):
    for i in range(N_ATOM - 1):
        # Вычисляется разность векторов между
        # i атомом и остальными, при этом i обрезается
        ijCoord = coords[i] - coords[i+1:]

        # Если расстояние между координатой атомов меньше половины ячейки
        # то расстояние увеличивается на ячейку
        # иначе уменьшается (берется переодический образ)
        ijCoord[ijCoord < -COUB_SIZE / 2] += COUB_SIZE
        ijCoord[ijCoord > COUB_SIZE / 2] -= COUB_SIZE

        # Вычисляем вспомог расстояния между атомами R
        dist2 = np.sum(ijCoord * ijCoord, axis=1)
        dist2Inv = 1.0 / dist2
        dist6Inv = dist2Inv * dist2Inv * dist2Inv

        enr = 4.0 * dist6Inv * (dist6Inv - 1.0)
        # вектор сил взаимодействия между i j атомами
        # ВЫчисляем вектор межатомных сил
        ff = 48.0 * dist2Inv * dist6Inv * (dist6Inv - 0.5)

        # Суммируем вектор межатомных сил с разностью координат
        # Получится n - i по 3
        # Каждое значение из ijCoord мы умножаем на строку в ff
        force_ = np.einsum('ki,k->ki', ijCoord, ff)

        # Добавляем i тую силу как сумму сил остальных
        forces[i] += np.sum(force_, axis=0)

        # Добавляем i в матрицу остальные силы
        forces[i+1:] -= force_
        
        energy += np.sum(enr)

    return [forces, energy]


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
scene.autoscale = False
gc1 = gcurve(color=color.cyan)
gc2 = gcurve(color=color.orange)


coords = initPositions()  # массив радиус-векторов частиц
acceleration = initAcceleration()  # массив ускорений
speed = initSpeed()  # массив скоростей
forces = initForces()  # массив векторов сил, действующих на каждую частицу

atoms = createAtomsByPos(coords)
energy = 0
temp = 0

# LIFECYCLE

for iter in range(ITERATIONS):
    #gc1.plot(pos=(iter, energy))
    gc2.plot(pos=(iter, temp))

    #scene.waitfor('click')

    [forces, energy] = calculation(coords, energy)
    temp = energy / 3.0 * N_ATOM * INTEGRATION_STEP ** 2
    
    # energy=energy+enr
    
    # вычисление (безразмерных) ускорений атомов системы
    acceleration = forces * INTEGRATION_STEP * INTEGRATION_STEP / 2
    
    # вычисление скоростей атомов
    speed += 2.0 * acceleration
    
    #print(coords[0],speed[0],acceleration[0])
    # новые положения частиц
    coords += speed + acceleration
    coords[coords < 0] += COUB_SIZE
    coords[coords > COUB_SIZE] -= COUB_SIZE

    
    
    changePosAtoms(atoms, coords)
# один раз за 5 (kstep) шагов результаты о координатах частиц выводятся в файл,
# и информация о номере шага выводится на экран
# для ускорения расчетов экзаменационного задания вывод в файл можно закомментировать

# Обнуляем массив сил
    forces[:, :] = 0
