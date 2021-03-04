import numpy as np
from typing import Tuple, List

# Linear Regresion
def linearReg(x, y):
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    x_sum = np.sum(x)
    y_sum = np.sum(y)
    x2_sum = np.sum(x*x)
    xy_sum = np.sum(x*y)

    m = ((n * xy_sum) - (x_sum * y_sum)) / ((n * x2_sum) - (x_sum * x_sum))
    b = (y_sum - (m * x_sum)) / n

    return m, b

def pbcPosition(s,L):
    if s > L:
        s -= L
    elif s < 0:
        s += L
    return s

def pbcSeparation(ds,L):
    if ds > L/2:
        ds -= L
    elif ds < -L/2:
        ds += L
    return ds

# Gets the center of mass of a finished system
# Formula found in article: Calculating Center of Mass in an Unbounded 2D Environment
    # Authors: Linge Bai and David E. Breen
    # DOI: 10.1080/2151237X.2008.10129266
def centerOfMass(L: int, particles: int, *args) -> tuple:
    if len(args) == 1:
        arr = args[0]
        x, y = np.hsplit(arr, 2)
        x = x.flatten()
        y = y.flatten()
    elif len(args) == 2:
        x = args[0]
        y = args[1]
    else:
        raise ValueError("Incorrect array entered for x and y.")

    epsilon_x, zeta_x, epsilon_y, zeta_y = 0, 0, 0, 0
    for (xi, yi) in zip(x, y):
        epsilon_x += np.cos((xi / L) * 2 * np.pi)
        zeta_x += np.sin((xi / L) * 2 * np.pi)
        epsilon_y += np.cos((yi / L) * 2 * np.pi)
        zeta_y += np.sin((yi / L) * 2 * np.pi)
    epsilon_x = epsilon_x / particles
    zeta_x = zeta_x / particles
    epsilon_y = epsilon_y / particles
    zeta_y = zeta_y / particles

    theta_x = np.arctan2(-zeta_x, -epsilon_x) + np.pi
    theta_y = np.arctan2(-zeta_y, -epsilon_y) + np.pi

    center_x = (theta_x / (2 * np.pi)) * L
    center_y = (theta_y / (2 * np.pi)) * L

    return center_x, center_y

#FIX
def moveToCenter(lat_size: int, x: List[float], y: List[float], cx = None, cy = None) -> Tuple[List, List]:
    if cx is None or cy is None:
        num_particles = len(x)
        cx, cy = centerOfMass(lat_size, num_particles, x, y)

    move_x = (lat_size / 2) - cx
    move_y = (lat_size / 2) - cy

    for i in range(len(x)):
        x[i] = x[i] + move_x
        y[i] = y[i] + move_y

        if x[i] >= lat_size:
            x[i] = x[i] - lat_size
        elif x[i] < 0:
            x[i] = lat_size + x[i]

        if y[i] >= lat_size:
            y[i] = y[i] - lat_size
        elif y[i] < 0:
            y[i] = lat_size + y[i]

    return x, y

def rg2(lat_size, particles, x, y):
    cx, cy = centerOfMass(lat_size, particles, x, y)
    print(cx, cy)
    rg_sum = 0
    for xi, yi in zip(x, y):
        dx = abs(xi - cx)
        dy = abs(yi - cy)
        if dx > (lat_size / 2):
            dx = lat_size - dx
        if dy > (lat_size / 2):
            dy = lat_size - dy
        rg_sum += (dx ** 2 + dy ** 2)

    Rg2 = (1 / particles) * rg_sum

    return Rg2

def rg2c(lat_size, mass1, mass2, cx1, cx2, rg1, rg2, cx_new, cy_new):
    dcx1, dcy1 = abs((cx_new) - cx1), abs((cy_new) - cy1)
    dcx2, dcy2 = abs((cx_new) - cx2), abs((cy_new) - cy2)

    if (dcx1 > (lat_size / 2)):
        dcx1 = lat_size - dcx1

    if (dcy1 > (lat_size / 2)):
        dcy1 = lat_size - dcy1


    if (dcx2 > (lat_size / 2)):
        dcx2 = lat_size - dcx2

    if (dcy2 > (lat_size / 2)):
        dcy2 = lat_size - dcy2

    h1 = pow(dcx1, 2) + pow(dcy1, 2)
    h2 = pow(dcx2, 2) + pow(dcy2, 2)

    rg2_cm11 = rg1 + h1
    rg2_cm12 = rg2 + h2

    rg2_f = ((mass1 * rg2_cm11) + (mass2 * rg2_cm12)) / (mass1 + mass2)

    return rg2_f


if __name__ == "__main__":
    pass