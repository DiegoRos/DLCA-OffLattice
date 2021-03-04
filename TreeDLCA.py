import sys
import numpy as np
import matplotlib.pyplot as plt
from FunctionsDLCA import moveToCenter

class DLCATree:
    def __init__(self, x: float, y: float):
        self.head_x = x
        self.head_y = y
        self.children = []

    def __str__(self, level=0):
        ret = f"{level}:" + "\t" * level + f"({self.head_x}, {self.head_y})" + "\n"
        level += 1
        for child in self.children:
            ret += child.__str__(level)

        return ret

    def addChild(self, tree_object):
        self.children.append(tree_object)

    def getPositions(self):
        x = self.getX()
        y = self.getY()

        return x, y

    def getX(self):
        x = [self.head_x]
        for child in self.children:
            x += child.getX()

        return x

    def getY(self):
        y = [self.head_y]
        for child in self.children:
            y += child.getY()

        return y


def createTree(lat_size: int, x: list, y: list, node: DLCATree, tolerance=4e-1):
    l_dist = lat_size - 1
    print("This many particles left: ", len(x))
    i = 0
    while i < len(x): # Iterate through positions
        dx = node.head_x - x[i]
        dy = node.head_y - y[i]

        dist = np.sqrt(np.square(dx) + np.square(dy))

        # Checks to see if fits requirements
        if (l_dist - tolerance <= dist <= l_dist + tolerance) and (((np.abs(dx) >= lat_size - 2) and (np.abs(dy) <= 1 + tolerance)) or ((np.abs(dy) >= lat_size - 2) and (np.abs(dx) <= 1 + tolerance))):
            # If requirements are true the x and y values are removed from the function
            node.addChild(DLCATree(x.pop(i), y.pop(i)))

        elif (1 - tolerance <= dist <= 1 + tolerance):
            node.addChild(DLCATree(x.pop(i), y.pop(i)))

        else:
            i += 1

    for child in node.children:
        createTree(lat_size, x, y, child)

def unfold(counter, lat_size: int, node: DLCATree, move_dist_x: int= 0, move_dist_y:int = 0):
    node.head_x += move_dist_x
    node.head_y += move_dist_y
    counter += 1
    print("Number of particles unfolded:", counter)

    for child in node.children:
        dx = node.head_x - child.head_x
        dy = node.head_y - child.head_y

        if (dx >= lat_size / 2) and (dy >= lat_size / 2):
            counter = unfold(counter, lat_size, child, move_dist_x=lat_size, move_dist_y=lat_size)
        elif (dx <= -lat_size / 2 and dy <= -lat_size / 2):
            counter = unfold(counter, lat_size, child, move_dist_x=-lat_size, move_dist_y=-lat_size)
        elif (dx <= -lat_size / 2 and dy >= lat_size / 2):
            counter = unfold(counter, lat_size, child, move_dist_x=-lat_size, move_dist_y=lat_size)
        elif (dx >= lat_size / 2 and dy <= -lat_size / 2):
            counter = unfold(counter, lat_size, child, move_dist_x=lat_size, move_dist_y=-lat_size)
        elif dx >= lat_size / 2:
            counter = unfold(counter, lat_size, child, move_dist_x=lat_size)
        elif dx <= -lat_size / 2:
            counter = unfold(counter, lat_size, child, move_dist_x=-lat_size)
        elif dy >= lat_size / 2:
            counter = unfold(counter, lat_size, child, move_dist_y=lat_size)
        elif dy <= -lat_size / 2:
            counter = unfold(counter, lat_size, child, move_dist_y=-lat_size)
        else:
            counter = unfold(counter, lat_size, child)

    return counter


def plot(x, y, lat_size, particles, original = False):
    x_diff = np.amax(x) - np.amin(x)
    y_diff = np.amax(y) - np.amin(y)

    if x_diff > y_diff:
        dist = x_diff + 10
    else:
        dist = y_diff + 10

    fig = plt.figure(figsize= (9,9))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title("DLCA Cluster L = {}, Particles = {}".format(lat_size, particles), fontsize = 16)
    rad = 1 / 2
    for xi, yi in zip(x, y):
        circ = plt.Circle((xi,yi), rad, color = '#002020')
        ax.add_patch(circ)
    if original:
        ax.set_xlim(0,lat_size)
        ax.set_ylim(0,lat_size)
    else:
        ax.set_xlim(np.min(x) - 10, np.min(x) + dist)
        ax.set_ylim(np.min(y) - 10, np.min(y) + dist)
    ax.set_xlabel("x", fontsize = 16)
    ax.set_ylabel("y", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)


if __name__ == "__main__":
    sys.setrecursionlimit(10000)

    lat_size = 50
    particles = 200
    cluster = np.loadtxt(f"Results/ClusterSize{lat_size}Particles{particles}.csv", delimiter=',')
    x, y, *_ = np.hsplit(cluster, cluster.shape[1])
    x = list(x.flatten())
    y = list(y.flatten())

    cluster_tree = DLCATree(x.pop(), y.pop())
    createTree(lat_size, x, y, cluster_tree)
    print("Tree Created \n")
    # print(cluster_tree)

    unfold(0, lat_size, cluster_tree)
    # print(cluster_tree)

    x1, y1 = cluster_tree.getPositions()
    cx = sum(x1) / len(x1)
    cy = sum(y1) / len(y1)
    Rg2 = (1 / particles) * sum([((xi - cx) ** 2) + ((yi - cy) ** 2) for xi, yi in zip(x1, y1)])
    print(Rg2)
    # print(x1)
    # print(y1)

    x, y, *_ = np.hsplit(cluster, cluster.shape[1])
    x = list(x.flatten())
    y = list(y.flatten())

    x, y = moveToCenter(lat_size, x, y)

    plot(x,y,lat_size,particles,original=True)
    plot(x1, y1, lat_size, particles)

    plt.show()