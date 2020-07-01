# МД-моделирование  Леннард-Джонсовой  жидкости

import time
import numpy as np
from vpython import *
scene = canvas(width=500, height=500)
scene.camera.pos = vector(3, 3, 3)
scene.autoscale = False

#  число атомов в системе (можно поменять на 2048, чтобы оценить проблемы с быстродействием)
N_ATOM = 256
# безразмерные: плотность, температура и шаг интегрирования системы по времени
DENSTY = 0.9
TEMP = 1.5
STEP = 0.001


def initArrays():
    # Раздел выделения памяти под основные массивы данных (без заполнения самими данными)
    # массив радиус-векторов частиц
    x_ = np.empty(N_ATOM*3).reshape(N_ATOM, 3)

    # массив ускорений
    a_ = np.empty(N_ATOM*3).reshape(N_ATOM, 3)

    # массив векторов сил, действующих на каждую частицу (инициализируется нулями)
    f_ = np.zeros(N_ATOM*3).reshape(N_ATOM, 3)
    return [x_, a_, f_]

# квадрат шага по времени и половина этой величины
# (используются при численном интегрировании)
stepsq = STEP*STEP
stepsqh = 0.5*stepsq

# расчет координат атомов основан на схеме их расположения
# в узлах гранецентрированной решетки
# (в экзаменационном задании схема расположения более простая! Плотность задавать не нужно!)

# Объем вычисляется на основе инфоромации о плотности и количестве атомов
vol = N_ATOM/DENSTY
# длина ребра куба - основной ячейки моделирования
cube = vol**(1/3)
cubeh = 0.5*cube

def initPositions(x_):
    # расчеты количества повторений блока решетки в пространстве для заполнения
    # объема моделирования

    # Код формально переведен с Фортрана.
    # Формат гранецентрированной решетки не допускает использование произвольного количества атомов:
    # их число только такое, которое может быть размещено в кубе. 256 и 2048 - подходят.
    # В экзаменационном задании эти расчеты не нужны!
    nunit = (N_ATOM/4.)**(1/3)+0.1
    ncheck = 4*(nunit**3)
    while ncheck < N_ATOM:
        nunit += 1
        ncheck = 4*(nunit**3)
    dist = 0.5*cube/nunit
    x_[0] = np.array([0.0, 0.0, 0.0])
    x_[1] = np.array([0.0, dist, dist])
    x_[2] = np.array([dist, 0.0, dist])
    x_[3] = np.array([dist, dist, 0.0])


    # Задание начальных положений атомов
    # В экзаменационной работе должна быть другая схема, в соответствии с вариантом А или Б
    m = 0
    kct = 0
    for i in range(int(nunit)):
        for j in range(int(nunit)):
            for k in range(int(nunit)):
                for ij in range(4):
                    if kct < N_ATOM:
                        x_[ij+m] = x_[ij]+2.0*dist * \
                            np.array([k, j, i], dtype=np.float64)
                        kct += 1
                m += 4
    return x_


x_, a_, f_ = initArrays()

x_ = initPositions(x_)

# Начальные скорости и их нормировка к начальной температуре
# равномерно распределенные случайные значения компонент скоростей частиц
v_ = np.random.uniform(-1, 1, (N_ATOM, 3))

# нормировка амплитуд скоростей, чтобы кинетическая энергия соотвествовала начальной температуре
velsq = np.sum(v_**2)
aheat = 3.e0*N_ATOM*stepsq*TEMP
factor = np.sqrt(aheat/velsq)
v_ *= factor




def createPoints(positions):
    points = []
    for p in positions:
        p = vector(p[0], p[1], p[2])
        points.append(sphere(pos=p, radius=0.1))
    return points

def changePosition(points, positions):
    for i in range(len(points)):
        p = positions[i]
        p = vector(p[0], p[1], p[2])
        points[i].pos = p
        
points = createPoints(x_)

# количество шагов по времени
# 1 кол-во шагов = 1 пентасекунда / шаг интегрирования 
maxeq = 1001
# итоговое время (step * maxeq) при переходе к размерным параметрам должно быть не меньше 1 пс
# в данном варианте - меньше!
# для подсчетов энергии
# в данном варианте программе не используется
energy = 0


# Цикл вычисления сил, действующих на частицу
start_time = time.time()
for ktime in range(maxeq):
    
    scene.waitfor('click')
    
    
    for i in range(N_ATOM-1):
        xij = x_[i] - x_[i+1:]
        xij[xij < -cubeh] += cube
        xij[xij > cubeh] -= cube
        rsq = np.sum(xij*xij, axis=1)
        rsqinv = 1.0/rsq
        r6inv = rsqinv*rsqinv*rsqinv
# формулу для вычисления потенциальной энергии необходимо скорректировать
# с учетом того, что r6inv в данной реализации - это вектор, а не скаляр, как было в Фортране
#        enr=4.0*r6inv*(r6inv-1.0)
        force_ = np.einsum('ki,k->ki', xij, rsqinv*48*r6inv*(r6inv-0.5))

        f_[i] += np.sum(force_, axis=0)
        f_[i+1:] -= force_
#        energy=energy+enr

# вычисление (безразмерных) ускорений атомов системы
    a_ = f_*stepsqh
    
# вычисление скоростей атомов
    v_ += 2.0*a_
# новые положения частиц
    x_ += v_+a_
    x_[x_ < 0] += cube
    x_[x_ > cube] -= cube

    changePosition(points, x_)
# один раз за 5 (kstep) шагов результаты о координатах частиц выводятся в файл,
# и информация о номере шага выводится на экран
# для ускорения расчетов экзаменационного задания вывод в файл можно закомментировать
    
# Обнуляем массив сил
    f_[:, :] = 0
# конец цикла по времени
delta_t = time.time() - start_time
print('dynamics successfully ended in', delta_t, 'seconds')
