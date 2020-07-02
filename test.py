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

FULL_TIME = 3.5e-5
INTEGRATION_STEP = 3.5e-8
ITERATIONS = int(FULL_TIME / INTEGRATION_STEP)


def initSpeed():
    """ Строит случайные начальные скорости """

    v_ = np.random.uniform(-1, 1, (N_ATOM, 3))
    # нормировка амплитуд скоростей, чтобы кинетическая энергия соотвествовала начальной температуре
    velsq = np.sum(v_ ** 2)
    
    aheat = 3.e0 * N_ATOM * INTEGRATION_STEP ** 2 * TEMP
    factor = np.sqrt(aheat / velsq)
    
    print(velsq)
    print(factor)

    return v_ * factor

print(initSpeed())