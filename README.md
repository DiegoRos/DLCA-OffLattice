# DLCA Off-Lattice Project
  ## Author: Diego Rosenberg

## Description 
This code generates a DLCA cluster through a random collocation of particles and a random movement of these particles (which is dependent on their mass) to generate larger clusters. The movement of these particles is pseodorandom and is dependant on their mass (as would happen in a particle dynamics simulation). Along with this a restriction can be added so that particles can only stick if they have a coordination number of 0 or 1. On these clusters we measure their fractal dimension and the probability of the cluster percolating through the lattice.

## Code
This project is divided into three main files in charge of creting the cluster, one for measuring the fractal dimension of resultant clusters and another file responsible for measuring the probability of the cluster percolating.

- Creation of CCA Cluster:
  - [x] MainDLCA.c
- Measurment of Fractal Dimension: 
  - [x] FracDimDLCA.py
- Measurment of Fractal Percolation:
  - [x] PercolatesDLCA.c
- Extras:
  - [x] Unfold Cluster &#8594; TreeDLCA.py
  - [x] Center of Mass and Move to Center in Python &#8594; FunctionsDLCA.py
  - [x] Plot Clusters in Python &#8594; DlCAPlot.py 

## Files
### MainCCA.c
Contains the main code and is responsible for generating the clusters.
- Imported Libraries:
  - <stdio.h>, <sys/stat.h>, <math.h>, <stdlib.h>, <time.h>
- Defined Constants:
  - Pi: 3.14159265358979
  - True: 1 
  - False: 0
  - Optional: MAX_COORDINATION (this will make a block of optional code be activated that will limit the corrdination number of all particles)
- Structures:
  - Particle
    - Variables:
      - double x, y &#8594; Particle position
      - int number, index &#8594; Particle number (I.D.) and index (cluster I.D.)
      - int coordination_number &#8594; Number of particles connected to it
      - int k &#8594; Cell position of particle
    - "Methods" (Not actually methods but sets the value of a variable inside an instance of the struct)
      - void setK
        - args: struct Particle *particle
        - Description: Sets k value of a particle struct using its x and y position.
  - Cluster
    - Variables:
      - int mass &#8594; Mass of cluster
      - float cx, cy &#8594; Center of mass of cluster
      - float rg2 &#8594; Radius of Gyration squared of cluster
- Functions:
  -  int main
    - args: int argc, char *argv[]
      - argv is necesary to execute the program in the command line (or as args if in a IDE), when running the program the folowing inputs are required as a minimum &#8594; &#60;Executable Name&#62; &#60;Number of Particles&#62; &#60;Lattice Size&#62;, while other two charactaristics can be altered &#60;Ring Size&#62; &#60;Progress&#62;, giving a total of 3 to 5 conditions that can be given to the program when running.
   - Description: This function will take in the inputs speficied in argv and create a completed (or in a state of equilibrium) DLCA cluster utilizing the following functions. This function prints the progress of the cluster to a file int the "Partial Results" directory and will print the finalized cluster to the "Results" directory. **(Note: Both of these directories are created at run time if they are not previously present)** The partial files have the name of PartialClusterSize& **#60;Lattice Size&#62;**Particles **&#60;Number of Particles&#62;**.csv while the final result has the form of ClusterSize **&#60;Lattice Size&#62;** Particles **&#60;Number of Particles&#62;**.csv
  - void allocate_memory
    - args: None
    - Description: Allocates memory to all of the pointer lists with the appropriate values.
  - void deallocate_memory
    - args: None
    - Description: Frees all memory created by the allocate_memory function.
  - void initialize
    - args: None
    - Description: First function ran in the procedure, creates all necesary structures and places all particles in the system (in such a way that none of the particles touch each other), as well as sets the mass_list, denominator and constant of normalization to initial values.
  - void setKlist
    - args: int k, int *kptr
    - Description: Function which returns the 9 cells that are neighbors of a given k cell. This is done so that the values are only calculated once throughout the entire procedure.
  - void addToCellList
    - args: struct Particle *particle
    - Description: Adds a particle's k to the Cell List to have a relative position of all the particles.
  - void resetCellListElement
    - args: struct Particle *particle
    - Description: Removes one particle from the cell list (to later be placed again after moving).
  - void changeDenominator
    - args: int mass1, int mass2
    - Description: When joining clusters the masses are joined to change the value of the denominator correspondingly.
  - void changeMassList
    - args: int mass1, int mass2
    - Description: When joining clusters the masses are joined and changed in the mass_list correspondingly.
  - void setA
    - args: None
    - Description: Sets the constant of normalization with the global denominator.
  - int checkSpot
    - args: double x, double y
    - Description: Function which takes a position and checks if this position is within interaction range of another particle. This function is utilized the first time particles are being placed in the system.
  - void step
    - args: int selected_cluster, double dir
    - Description: Given a cluster and a direction the cluster is moved in a specified direction, as well as changing cell_list and center of mass for the cluster accordingly.
  - void stepBack
    - args: int selected_cluster, double *odist, double dir
    - Description: Given a cluster, a direction, and a return distance the cluster is moved in a specified direction to remove overlap, as well as changing cell_list and center of mass of the cluster accordingly.
  - double periodicBoundaryConditions
    - args: double val
    - Description: If a particle is placed outside the periodic boundary conditions the function pushes it back on the opposite side.
  - int checkCluster
    - args: int selected_cluster, double dir, int *clusters, double *odist
    - Description: Given a cluster all values from *cluster and *odist are set to base values to then be compared, checks for every particle in the cluster.
  - void checkParticle
    - args: int selected_particle, double dir, int *checked_clusters, double *odist
    - Description: Chekcs if a particle from another cluster is in interaction range from it
  - float overlapDist
    - args: double dist, double dir, double p1x, double p1y, double p2x, double p2y
    - Description: Calculates the overlap distance of two neigboring particles. The formula is created utilizing a composition of triangles and other geometric tricks.
  - void joinClusters
    - args: int c1, int c2
    - Description: Function to join two clusters by index
  - void posParticlesCluster
    - args: int lc_label, double *x_list, double *y_list
    - Description: Function that gives the position of all particles in a cluster utilizing the linked list
  - void centerOfMass
    - args: int particles, double *x_list, double *y_list, double *cx, double *cy
    - Description: Calculates Center of mass for a given x_list and y_list in periodic boundary conditions. Formula found in article: Calculating Center of Mass in an Unbounded 2D Environment, Authors: Linge Bai and David E. Breen, DOI: 10.1080/2151237X.2008.10129266

### FracDimDLCA.c
- This code will calculate the fractal dimension for a 2D cluster that is off-lattice or on-lattice. This will be done utilizing two methods, from these the average will represent the fractal dimension. The methods are:
  1. The absolute value of the slope of the linear regression between logarithm of the box size vs logarithm the average particles per box.
  2. The absolute value of the slope of the linear regression of the logarithm box size vs logarithm number of boxes.
- Imoprted Libraries:
  - &#60;stdio.h&#62;, &#60;sys/stat.h&#62;, &#60;math.h&#62;, &#60;stdlib.h&#62;, &#60;time.h&#62;, &#60;string.h&#62;
- Functions:
  - int main
    - args: int argc, char *argv[]
      - When executing the compiled file there have to be three entries without including the executable name, resulting in 4 total entries, these are: &#60;Executable Name&#62; &#60;Lattice Size&#62; &#60;Number of Particles&#62; &#60;File Name&#62;. 
    - Description: Taking the argv entries the size of the box_counts, avg_counts and box_size pointer lists are calculated, then the fractal dimension function function is called and the returned value, along with the resulting lists (not the logaritmic values) are saved to a file in the Results directory with the format: FracDimDCountsSize **&#60;Lattice Size&#62;** Particles **&#60;Number of Particles&#62;**.csv and the fractal dimension is saved to: FracDimSize **&#60;Lattice Size&#62;** Particles **&#60;Number of Particles&#62;**.txt
  - void allocate_memory
    - args: None
    - Description: Allocates memory to all of the pointer lists with the appropriate values.
  - void deallocate_memory
    - args: None
    - Description: Frees all memory created by the allocate_memory function.
  - double fractalDimension
    - args: double *x, double *y
    - Description: Calculates the fractal dimension utilizing the global pointer lists named box_counts, avg_counts and box_size.
  - void fractalDimensionCounts
    - args: double *x, double *y
    - Perfomrs the box count and average count methods to the lists of x and y sent in, the results are stored in the global pointer lists box_counts, avg_counts and box_size.
  - void linearReg
    - args: int n, double *x, double *y, double *m, double *b
    - Description: returns 2 floating point values which are the slope and intercept respectively, these values are returned in the m and b pointers.

### PercolatesDLCA.c
- This code will check every row and every column of the system to revise if a given cluster percolates the system, if all of the columns are filled by at least one particle the cluster will percolate in x, and if all rows are filled by at least one particle the cluster will percolate in the y direction.

### FunctionsCCA.py
Contains the extra functions necessary to find the center of mass and center a given cluster.
- Imported Libraries:
  - numpy, typing
- Functions: 
    - centerOfMass
      - **args**: L, particles, *args (list of x and y values or numpy array of x and y values)
      - **kwargs**: NONE
      - This function calculates the center of mass for a cluster existing inside a lattice with periodic boudary conditions utilizing the formula found in article:Calculating Center of Mass in an Unbounded 2D Environment. Authors: Linge Bai and David E. Breen. DOI: 10.1080/2151237X.2008.10129266
      - returns two floating values representing the center of mass in x and the center of mass in y.
    - moveToCenter
      - **args**: lat_size, x, y
      - **kwargs**: cx (Default value is None), cy (Default value is None)
      - Moves the center of mass of the cluster to the center of the lattice, this can be usefull to plot clusters in a better way.

