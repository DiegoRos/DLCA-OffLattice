"""
========================================================================================================================
Code: Multiple DLCA Cluster generator in 2D
File: MainDLCA.py
Author: Diego Rosenberg

Description: This code generates multiple DLCA clusters through utilizing MainDLCA.py, this saves the resulting
clusters in a given directory and saves the results of the time elapsed, number of particles, lattice size, fractal
dimension, if the cluster percolates or not, as well as all of the coordination number for the entire cluster.
========================================================================================================================
"""
import os
import sys
import subprocess
from shutil import copyfile
import datetime
import time
import csv
import pandas as pd
import numpy as np


def mainScript(lat_size, particles, total, make=True):
    result_file = f"Results/ClusterSize{lat_size}Particles{particles}.csv"

    if make == True:
        os.system("make")

    main_executable_name = "DLCA "
    fracdim_executable_name = "FracDimDLCA "
    percolation_executable_name = "Percolates "

    field_names = ["date_and_time", "elapsed_time_sec", "lattice_size", "num_particles", "fractal_dim", "rg2","percolates",
                   "z_1", "z_2", "z_3", "z_4", "z_5", "z_6"]
    # Create results directory to store a .csv containing results from all clusters.
    if not os.path.isdir('Final Results'):
        os.mkdir('Final Results')

    # Creates .csv file for specific cluster lattice size and amount of particles, will store date run, elapsed time,
    # lattice size, number of particles, fractal dimension and if the cluster percolates the system.
    if not os.path.isfile('Final Results/ResultsOffL.csv'):
        with open('Final Results/ResultsOffL.csv', mode='w') as file:
            writer = csv.DictWriter(file, fieldnames=field_names)
            writer.writeheader()

    # If a previous file existed this will count the number of lines to only add the amount of required rows.
    df = pd.read_csv('Final Results/ResultsOffL.csv')
    try:
        num_lines = len(df[df['num_particles'] == particles])
    except:
        num_lines = 0

    results_sims = []
    # Runs the main loop and fractal dimension code n times and writes it to .csv to keep results.
    for i in range(num_lines, total):
        with open('Final Results/ResultsOffL.csv', mode='+a', newline='') as file:
            dict_writer = csv.DictWriter(file, fieldnames=field_names)

            date = datetime.datetime.now()  # Gets current date and time
            start = time.time()  # Starts internal clock

            # Run cluster executable:
            os.system(main_executable_name + f"{lat_size} {particles} 10000")

            cluster = np.loadtxt(result_file, delimiter=',')
            x, y, indexes, coordination_numbers = np.hsplit(cluster, cluster.shape[1])
            x = x.flatten()
            y = y.flatten()

            with open(f"Results/Rg2ClusterSize{lat_size}Particles{particles}.txt", 'r') as frg:
                rg2 = float(frg.readline())

            coordination_numbers = coordination_numbers.flatten()
            *_, coordination_count = np.unique(np.array(coordination_numbers), return_counts=True)
            for _ in range(6 - len(coordination_count)):
                coordination_count = np.append(coordination_count, 0)

            percolation = bool(os.system(percolation_executable_name + f"{lat_size} {particles} \"{result_file}\""))

            # Calculates fractal dimension
            os.system(fracdim_executable_name + f"{lat_size} {particles} \"{result_file}\"")
            with open(f"Results/FracDimSize{lat_size}Particles{particles}.txt", 'r') as fd:
                fractal_dimension = float(fd.readline())

            end = time.time()  # End time for internal clock

            # Writes all values to data base
            print('Writing')
            dict_of_elements = dict(zip(field_names, [
                                    str(date), (end - start), lat_size, particles, fractal_dimension, rg2, percolation,
                                    int(coordination_count[0]), int(coordination_count[1]), int(coordination_count[2]), 
                                    int(coordination_count[3]), int(coordination_count[4]), int(coordination_count[5])
                                    ]
                                ))

            dict_writer.writerow(dict_of_elements)
            print(dict_of_elements)

            np.savetxt(f"D:\\ASE III Resultados/{i}ClusterSize{lat_size}Particles{particles}.csv",
                        np.column_stack([x, y, coordination_numbers]), delimiter=',')

            fd_lists = np.loadtxt(f"Results/FracDimCountsSize{lat_size}Particles{particles}.csv", delimiter=',')
            box_sizes, box_counts, avg_counts = np.hsplit(fd_lists, fd_lists.shape[1])
            box_sizes = box_sizes.flatten()
            box_counts = box_counts.flatten()
            avg_counts = avg_counts.flatten()

            np.savetxt(f"D:\\ASE III Resultados/{i}FracDimListSize{lat_size}Particles{particles}.csv",
                        np.column_stack([box_sizes, box_counts, avg_counts]), delimiter=',')

            copyfile(f"Results/RgMassTimeSize{lat_size}Particles{particles}.csv", f"D:\\ASE III Resultados/{i}RgMassTimeSize{lat_size}Particles{particles}.csv")
            

if __name__=='__main__':
    if len(sys.argv) == 4:
        lat_size = int(sys.argv[1])
        particles = int(sys.argv[2])
        total = int(sys.argv[3])
        mainScript(lat_size, particles, total)

    else:
        print("Incorrect number of entries, read documentation.")
