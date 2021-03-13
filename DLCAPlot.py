import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from FunctionsDLCA import *

def plot(lat_size, particles, file_name, center = False):
    cluster = np.loadtxt(file_name, delimiter=',')

    x, y, indexes, *_ = np.hsplit(cluster, cluster.shape[1])
    x = x.flatten()
    y = y.flatten()
    indexes = indexes.flatten()

    if center:
        x, y = moveToCenter(lat_size, x, y)

    cividis = cm.get_cmap("cividis", int(np.amax(indexes))+1)

    fig = plt.figure(figsize= (9,9))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title("DLCA Cluster L = {}, Particles = {}".format(lat_size, particles), fontsize = 16)
    rad = 1 / 2
    for xi, yi, index in zip(x, y,indexes):
        circ = plt.Circle((xi,yi), rad, color = cividis(index / (np.amax(indexes) + 1)))
        ax.add_patch(circ)
    ax.set_xlim(0, lat_size)
    ax.set_ylim(0, lat_size)
    ax.set_xlabel("x", fontsize = 16)
    ax.set_ylabel("y", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        lat_size = int(sys.argv[1])
        particles = int(sys.argv[2])
        chosen_file = '-m'

    elif (len(sys.argv) == 4):
        lat_size = int(sys.argv[1])
        particles = int(sys.argv[2])
        chosen_file = str(sys.argv[3])

    else:
        quit()

    if chosen_file == '-m':
        file_name = f"Results/ClusterSize{lat_size}Particles{particles}.csv"

    elif chosen_file == '-p':
        file_name = f"Partial Results/PartialClusterSize{lat_size}Particles{particles}.csv"

    elif chosen_file == '-me':
        file_name = f"Results/EdgeClusterSize{lat_size}Particles{particles}.csv"

    elif chosen_file == '-pe':
        file_name = f"Partial Results/EdgePartialClusterSize{lat_size}Particles{particles}.csv"
    
    else:
        quit()

    plot(lat_size, particles, file_name, center=True)

    plt.show()
