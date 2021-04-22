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

# Function that counts the average number of particles per box in different divisions, calculates the size of these boxes,
# and calculates the amount of boxes that contain particles for a given division size.
def boxCount(L: int, particles: int, x: List[float], y: List[float]) -> Tuple[list, list, list]:
    div = 4 + (4 * (particles <= 8000)) # Sets start division
    break_condition = 2 + (2 * (L > 100)) # Sets ending place for amount of divisions
    box_count = [] # List containing the amount of boxes with particles
    avg_counts = [] # List containing the amount of particles average per box
    box_size = [] # List containing the size of the current boxes
    while True: # Emulation of do while loop to calculate individual vales to put into previous lists
        boxcount = 0
        division = L / div

        if division < break_condition: # Ends do while loop
            break

        for i in range(1,div+1):
            for j in range(1,div+1):
                lower_x = ((i - 1) * division)
                upper_x = (i * division)
                lower_y = ((j-1) * division)
                upper_y = (j * division)
                for xi, yi in zip(x, y):
                    if lower_x <= xi <= upper_x and lower_y <= yi <= upper_y:
                        boxcount += 1
                        break # If a box contains 1 particle the code can move on to the next

        # Appends found values to the lists
        box_count.append(boxcount)
        avg_counts.append(particles/boxcount)
        box_size.append(division)

        print("{:.3f} divisions with minimum {}".format(division, break_condition))

        div *= 2

    return box_count, avg_counts, box_size


if __name__ == "__main__":
    pass