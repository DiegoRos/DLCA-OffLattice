import os
import sys
import subprocess
import shutil
import re
import numpy as np
import imageio

import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams.update({'figure.max_open_warning': 0})
plt.rcParams.update({'font.size': 16})

global cluster_fn
cluster_fn = "Cluster"
global animation_folder
animation_folder = "Animation/"

def mainGif(lat_size, partilces, probability=0, run_sim=True):

    if run_sim:
        val = os.system("gcc -D GIF DLCA.c -o DLCAGif.exe -lm")
        if val != 0:
            cont = print("Could not compile DLCA.c file with Gif.")
            exit()

        try:
            shutil.rmtree('Animation') # Eliminates existing path with (possible) previous images, will not work if folder open
        except:
            print()

        progress = int(input("Enter progress values: "))
        if progress:
            subprocess.run(["DLCAGif.exe", str(lat_size), str(particles), str(probability), str(progress)])
        else:
            subprocess.run(["DLCAGif.exe", str(lat_size), str(particles), str(probability)])

    print("Creating Images")

    text_files = [f for f in os.listdir("Animation") if (cluster_fn in f) and (str(int(lat_size)) in f) and (str(int(particles)) in f) and (str(float(probability)) in f) and (f.endswith(".csv"))]

    text_files = sorted(text_files, key=lambda x: int(x.split('-')[0]))
    for i, file_i in enumerate(text_files):

        cluster = np.loadtxt(animation_folder + file_i, delimiter = ',')

        x, y, index = np.hsplit(cluster, cluster.shape[1])
        x = x.flatten()
        y = y.flatten()
        indexes = index.flatten().astype(int)
        max_indexes = int(np.max(indexes)) + 1

        cividis = cm.get_cmap("cividis", max_indexes)
        fig = plt.figure(figsize=(9, 9))
        ax = fig.add_subplot(1, 1, 1)
        if probability != 0:
            ax.set_title(f"DLCA Reversible Cluster Lattice Size = {lat_size}, \nParticles = {particles} & Probability = {probability}")
        else:
            ax.set_title(f"DLCA Cluster Lattice Size = {lat_size}, Particles = {particles}")
        rad = 1 / 2
        
        for xi, yi, index in zip(x, y, indexes):
            circ = plt.Circle((xi, yi), rad, color = cividis(index / max_indexes))
            ax.add_patch(circ)
        ax.set_xlim(0, lat_size)
        ax.set_ylim(0, lat_size)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.savefig(f"Animation/{i}-ClusterSize{lat_size}Particles{particles}Prob{probability}.png", dpi=200)
        plt.close()

        os.remove(animation_folder + file_i)

    print("Creating Movie")
    with imageio.get_writer(f'Animation/MovieSize{lat_size}Particles{particles}Prob{probability}.gif', mode='I') as writer:
        if probability != 0:
            pictures = [f for f in os.listdir("Animation") if (cluster_fn in f) and (str(int(lat_size)) in f) and (str(int(particles)) in f) and (str(float(probability)) in f) and (f.endswith(".png"))]
        else:
            pictures = [f for f in os.listdir("Animation") if (cluster_fn in f) and (str(int(lat_size)) in f) and (str(int(particles)) in f) and (str(int(probability)) in f) and (f.endswith(".png"))]

        pictures = sorted(pictures, key=lambda x: int(x.split('-')[0]))
        for filename in pictures[:-1]:
            image = imageio.imread(animation_folder + filename)
            writer.append_data(image)
            os.remove(animation_folder + filename)

        image = imageio.imread(animation_folder + pictures[-1])
        for i in range(15):
            writer.append_data(image)
        os.remove(animation_folder + pictures[-1])



if __name__ == '__main__':
    if (len(sys.argv) < 3):
        exit()

    lat_size = int(sys.argv[1])
    particles = int(sys.argv[2])


    if '-p' in sys.argv:
        index = sys.argv.index('-p')
        probability = sys.argv[index + 1]
    else:
        probability = 0
    
    if '-r' in sys.argv:
        index = sys.argv.index('-r')
        run_sim = True
    else:
        run_sim = False

    mainGif(lat_size, particles, probability=probability, run_sim=run_sim)

