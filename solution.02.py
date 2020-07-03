# Периодические граничные условия (бесконечная среда)

# IMPORTS
import time
import numpy as np

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

FULL_TIME = 3e-5
INTEGRATION_STEP = 3e-8
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
    """ Строит нулевой массив ускорений """
    return np.zeros(N_ATOM*3).reshape(N_ATOM, 3)


def initSpeed():
    """ Строит случайные начальные скорости """
    v_ = np.random.uniform(-1, 1, (N_ATOM, 3))
    
    # нормировка амплитуд скоростей, чтобы кинетическая энергия соотвествовала начальной температуре
    velsq = np.sum(v_ ** 2)
    aheat = 3.0 * N_ATOM * INTEGRATION_STEP ** 2 * TEMP
    factor = np.sqrt(aheat / velsq)

    return  v_ * factor


def initForces():
    """ Строит нулевой массив векторов сил """
    return np.zeros(N_ATOM*3).reshape(N_ATOM, 3)


def evaluateForces(coords, forces, energy = 0):
    """ Вычисление сил действующих на частицы | вычисление энергии"""
    
    for i in range(N_ATOM - 1):
        # Вычисляется разность между радиус векторами частиц
        ijCoord = coords[i] - coords[i+1:]        

        # Вычисляем вспомог расстояния между атомами
        dist2 = np.sum(ijCoord * ijCoord, axis=1)
        dist2Inv = 1.0 / dist2
        dist6Inv = dist2Inv * dist2Inv * dist2Inv
        
        # Вычисляем межатомные силы
        ff = 48.0 * dist2Inv * dist6Inv * (dist6Inv - 0.5)
        force_ = np.einsum('ki,k->ki', ijCoord, ff)
        forces[i] += np.sum(force_, axis=0)
        forces[i+1:] -= force_
        
        # Потенциальная энергия Леннарда-Джонса
        enr = 4.0 * dist6Inv * (dist6Inv - 1.0)
        energy += np.sum(enr) 

    return energy

def werleScheme(acceleration, speed, coords, forces):
    """ Решение уравнений движения по 3й схеме Верле """
    # acceleration = forces * INTEGRATION_STEP * INTEGRATION_STEP / 2
    # speed += 2.0 * acceleration
    
    speed += (acceleration + forces) / 2 * INTEGRATION_STEP
    acceleration = forces
    forces[:, :] = 0
    
    coords += speed * INTEGRATION_STEP + acceleration * INTEGRATION_STEP * INTEGRATION_STEP / 2


# PREPARATION
dataCoordinates = open('data//data.coordinates.02.txt', mode='w')
dataPlot = open('data//data.plot.02.txt', mode='w')
print("Total energy of the system", file=dataPlot)

coords = initPositions()  # массив радиус-векторов частиц
acceleration = initAcceleration()  # массив ускорений
speed = initSpeed()  # массив скоростей
forces = initForces()  # массив векторов сил, действующих на каждую частицу

# LIFECYCLE
start_time = time.time()
for iter in range(ITERATIONS):
    
    """ Вычисление сил действующих на частицы + потенциальная энергия """
    energy = evaluateForces(coords, forces) 
    
    """ Решение уравнений движения / 3 схема Верле"""
    werleScheme(acceleration, speed, coords, forces)
    
    # Отталкивание c изменением температуры
    for i in range(N_ATOM - 1):
        for j in range(3):
            if (coords[i][j] <= 0) or (coords[i][j] >= COUB_SIZE):
                speed[i] *= -np.sqrt(TEMP_0 / TEMP_0)

    energy += np.sum(speed ** 2)
    
    if iter % 50 == 0: print(iter)
    
    for c in coords: print(c[0], c[1], c[2], file=dataCoordinates)
    print('', file=dataCoordinates)
    print(energy, file=dataPlot)

delta_t = time.time() - start_time
print('dynamics successfully ended in', delta_t, 'seconds') # 375 seconds

dataCoordinates.close()
dataPlot.close()