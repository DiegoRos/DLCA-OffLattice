// ============================================================
// Code: DLCA 2D Cluster Generator
// File: DLCA.c
// Author: Diego Rosenberg
/*
Description: This code generates a 2D DLCA cluster with no 
restrictions. These particles have a diffusion coefficient and
move randomly withi the system.
*/
// ============================================================

#include "DLCA.h"

int num_particles; // Amount of particles in the system
int lat_size; // Size of lattice of the system
int ring_size; // Diameter of particles and the distance of interaction between them
int progress; // Gives the amount of steps required to save a progress csv and to print an update to the console
float prob_separate;
int steps_taken = 0;
int break_condition;


int number_of_clusters; // The current number of clusters in the system

int *firstp; // Size of number of clusters
int *nextp; // Linked list with size of number of particles
int *lastp; // Size of number of cluster
Cluster *cluster_list; // Size of number of clusters

Particle *particle_list; // Size of number of particles

int cells; // Amount of 
int cells2;
int *cell_list; // Size of cells2 + number of particles

#define k_next 9
int **k_list;

int *mass_list; // Size of number of particles
double denominator;
double A;

Z1 list_z1;

// Instance methods for structs:

    // "Instance methods" for struct Particle (functions applied as if they were a instance method from OOP).
void setK(Particle *particle){
    (*particle).k = ((int) ((*particle).y / ring_size) * cells) + ((int) (*particle).x / ring_size); // Sets k value of a particle struct
}

// End of "instance methods"

// Program Function Declaration

void runSim();
void saveResults();
void allocate_memory();
void deallocate_memory();
void initialize();
void setKlist(int k, int *kptr);
void addToCellList(Particle *particle);
void resetCellListElement(Particle *particle);
void resetParticleInLists(int number);
void removeParticleZ1Lists(int number);
void addParticleZ1List(int number);
void changeDenominator(int mass1, int mass2);
void changeMassList(int mass1, int mass2);
void setA();
int selectClusterMass();
int checkSpot(double x, double y);
void step(int selected_cluster, double dir);
void stepBack(int selected_cluster, double *odist, double dir);
double periodicBoundaryConditions(double val);
int checkCluster(int selected_cluster, double dir, int *clusters, double *odist);
void checkParticle(int selected_particle, double dir, int *checked_clusters, double *odist);
float overlapDist(double dist, double dir, double p1x, double p1y, double p2x, double p2y);
void joinClusters(int c1, int c2);
int separateCluster(int number);
void posParticlesCluster(int lc_label, double *x_list, double *y_list);
void centerOfMass(int particles, double *x_list, double *y_list, int *mass_list, double *cx, double *cy);
void centerOfMassTwoPoints(double *cx_list, double *cy_list, int *mass_list, double *cx, double *cy);
void undoCenterOfMassTwoPoints(double *cx_list, double *cy_list, int *mass_list, double *cx, double *cy);
double radiusOfGyration(int lc_mass, int sc_mass, int lc_label, int sc_label, double *cx_new, double *cy_new);
double undoRadiusOfGyration(int sc_mass, int lc_label, int removed_particle, double *cx_new, double *cy_new);
void resetRgFile();
void writeRgFile(int cluster_index);
void resetNumberFile();
void writeNumberFile();

int main(int argc, char *argv[]){
    double time_spent = 0.0;
 
    clock_t begin = clock();

    // Options for argv when running executable

    if (argc < 3){
        printf("Not enough arguments entered when running DLCA.c executable.\n");
        printf("When running executable in cmd the order of entries is as follows: \n");
        printf("<Executable Name> <number of particles> <lattice size> <ring size> <progress>\n");
        exit(0);
    }

    else if (argc == 3){
        lat_size = atoi(argv[1]);
        num_particles = atoi(argv[2]);
        prob_separate = 0;
        progress = num_particles > 20000 ? 1000000 : (num_particles > 5000 ? 100000 : (num_particles));

		if (progress < 10){
			progress = 10;
		}

        ring_size = 1;
    }

    else if (argc == 4){
        lat_size = atoi(argv[1]);
        num_particles = atoi(argv[2]);
        prob_separate = atof(argv[3]);
        progress = num_particles > 2000 ? 1000000 : (num_particles > 500 ? 100000 : (num_particles));

        if (progress < 10){
            progress = 10;
        }

        ring_size = 1;
    }

    else if (argc == 5){
        lat_size = atoi(argv[1]);
        num_particles = atoi(argv[2]);
        prob_separate = atof(argv[3]);
        progress = atoi(argv[4]);
        ring_size = 1;
    }

    else if (argc == 6){
        lat_size = atoi(argv[1]);
        num_particles = atoi(argv[2]);
        prob_separate = atof(argv[3]);
        progress = atoi(argv[4]);
        ring_size = atoi(argv[5]);
    }

    else{
        printf("Too many arguments entered when running DLCA.c executable.\n");
        printf("When running executable in cmd the order of entries is as follows: \n");
        printf("<Executable Name> <number of particles> <lattice size> <probability> <progress> <ring size> \n");
        exit(0);
    }


    #ifdef CLUSTER_NUMBER
        resetNumberFile();
    #endif

    #ifdef DEBUG
        srand(0);
	#else
    	srand((unsigned)time(NULL)); // Set random seed 
    #endif


    printf("Allocating Memory...\n");
    allocate_memory();

    printf("Initializing Clusters...\n");
    initialize();

    printf("Running Simulation...\n");
    runSim();

    // Save results
    saveResults();
    
    deallocate_memory(); // Deallocate memory of global pointer variables

    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    printf("Time elpased is %f seconds\n", time_spent);

    return number_of_clusters;
}

void runSim(){
    // Create partial result directory if directory does not exist
    struct stat pr = {0};
    if (stat("/Partial Results", &pr) == -1){
        #ifdef __unix__
            mkdir("Partial Results", 0777);
        #else
            mkdir("Partial Results");
        #endif
    }
    char middle_file_name[80];
    #ifdef MAX_COORDINATION
        sprintf(middle_file_name, "Partial Results/EdgePartialClusterSize%dParticles%dProb%f.csv", lat_size, num_particles,prob_separate);
    #else
        sprintf(middle_file_name, "Partial Results/PartialClusterSize%dParticles%dProb%f.csv", lat_size, num_particles,prob_separate);
    #endif

    // Creation of animation files if required
    #ifdef GIF
        char gif_file_name[80];
        FILE *gif;
        // Create animation directory if directory does not exist
            // This will contain all partial results advancements in csv format to then be processed and graphed by python.
        struct stat a = {0};
        if (stat("/Animation", &a) == -1){
            #ifdef __unix__
                mkdir("Animation", 0777);
            #else
                mkdir("Animation");
            #endif
        }
    #endif

    FILE *fm;

    int selected_cluster, separation = False, cluster_one_counter = 0; // Randomly selected cluster index and boolean to revise if cluster separated
    int rand_z1;
    double rand_val;
    double dir; // Direction of step, will be changed every time step is called
    int *connecting_clusters = (int *)malloc(sizeof(int) * 4);
    int connected; // Boolean that checks if cluster connected
    double *odist = (double *)malloc(sizeof(double));
    steps_taken = 0;
    break_condition = ((prob_separate <= 0.15) && num_particles < 10000) ? (int)(lat_size * 1000000) : (int) (5 * lat_size * 1000000);

    //while(number_of_clusters != 1){
    while(True){ // Will run until break condition is reached
        if((prob_separate != 0) && (list_z1.total != 0)){//Condition to speed up code if the probability is 0
            rand_val = (double)rand() / (double)RAND_MAX; // Roll dice

            if(rand_val < prob_separate){
                separation = True;
                rand_z1 = rand() % (list_z1.total);

                int vecino = separateCluster(get(list_z1.start, rand_z1));
                selected_cluster = number_of_clusters - 1;
                double dx = particle_list[firstp[selected_cluster]].x - particle_list[vecino].x;
                double dy = particle_list[firstp[selected_cluster]].y - particle_list[vecino].y;
                if (dx > (lat_size / 2))
                    dx -= lat_size;
                if (dx < (-lat_size / 2))
                    dx += lat_size;

                if (dy > (lat_size / 2))
                    dy -= lat_size;
                if (dy < (-lat_size / 2))
                    dy += lat_size;
                dir = atan2(dy, dx); // Random value from 0 to 2 Pi that indicates the direction of movement.
            }
        }

        if (!separation){
            selected_cluster = selectClusterMass();
            dir = ((double)rand() / (double)RAND_MAX) * 2 * Pi; // Random value from 0 to 2 Pi that indicates the direction of movement.   
        }

        step(selected_cluster, dir);

        if (number_of_clusters != 1)
            connected = checkCluster(selected_cluster, dir, connecting_clusters, odist);

        #ifdef MAX_COORDINATION
            if(connected){
                if ((particle_list[connecting_clusters[0]].coordination_number < MAX_COORDINATION) && (particle_list[connecting_clusters[2]].coordination_number < MAX_COORDINATION)){
                   // Before joining clusters are pushed back and the overlap is removed
                    stepBack(selected_cluster, odist, dir);
                    if(prob_separate != 0){
                        if(particle_list[connecting_clusters[0]].coordination_number == 1){
                            removeParticleZ1Lists(connecting_clusters[0]);
                        }

                        if(particle_list[connecting_clusters[2]].coordination_number == 1){
                            removeParticleZ1Lists(connecting_clusters[2]);
                        }
                    }

                    // Add 1 to the coordination particle of the corresponding particles
                    particle_list[connecting_clusters[0]].coordination_number++;
                    particle_list[connecting_clusters[2]].coordination_number++;

                    if(prob_separate != 0){
                        if(particle_list[connecting_clusters[0]].coordination_number == 1){
                            addParticleZ1List(connecting_clusters[0]);
                        } 

                        if(particle_list[connecting_clusters[2]].coordination_number == 1){
                            addParticleZ1List(connecting_clusters[2]);
                        }
                    }

                    particle_list[connecting_clusters[0]].neighbor = push(connecting_clusters[2], particle_list[connecting_clusters[0]].neighbor);
                    particle_list[connecting_clusters[2]].neighbor = push(connecting_clusters[0], particle_list[connecting_clusters[2]].neighbor);

                    // join non overlaping cluster;
                    joinClusters(connecting_clusters[1], connecting_clusters[3]); 
                }

                else{
                    *odist = 0;
                    stepBack(selected_cluster, odist, dir);
                    steps_taken--;
                }
            }

        #else
            if(connected){
                // Before joining clusters are pushed back and the overlap is removed
                stepBack(selected_cluster, odist, dir);
                if(prob_separate != 0){
                    if(particle_list[connecting_clusters[0]].coordination_number == 1){
                        removeParticleZ1Lists(connecting_clusters[0]);
                    }

                    if(particle_list[connecting_clusters[2]].coordination_number == 1){
                        removeParticleZ1Lists(connecting_clusters[2]);
                    }
                }

                // Add 1 to the coordination particle of the corresponding particles
                particle_list[connecting_clusters[0]].coordination_number++;
                particle_list[connecting_clusters[2]].coordination_number++;

                if(prob_separate != 0){
                    if(particle_list[connecting_clusters[0]].coordination_number == 1){
                        addParticleZ1List(connecting_clusters[0]);
                    } 

                    if(particle_list[connecting_clusters[2]].coordination_number == 1){
                        addParticleZ1List(connecting_clusters[2]);
                    }
                }

                particle_list[connecting_clusters[0]].neighbor = push(connecting_clusters[2], particle_list[connecting_clusters[0]].neighbor);
                particle_list[connecting_clusters[2]].neighbor = push(connecting_clusters[0], particle_list[connecting_clusters[2]].neighbor);

                // join non overlaping cluster;
                joinClusters(connecting_clusters[1], connecting_clusters[3]);
            }

        #endif 

        if ((steps_taken % progress) == 0){
        // Print and save progress of system
            printf("Total number of clusters %d\n", number_of_clusters);

            fm = fopen(middle_file_name, "w");

            for(int i = 0; i < num_particles; i++){
                fprintf(fm, "%lf,%lf,%d,%d\n", particle_list[i].x, particle_list[i].y, particle_list[i].index, particle_list[i].coordination_number);
            }

            fclose(fm);

            #ifdef CLUSTER_NUMBER
                writeNumberFile();
            #endif

            #ifdef GIF
                sprintf(gif_file_name, "Animation/%d-ClusterSize%dParticles%dProb%f.csv", steps_taken, lat_size, num_particles, prob_separate);
                gif = fopen(gif_file_name, "w");

                for(int i = 0; i < num_particles; i++){
                    fprintf(gif, "%lf,%lf,%d\n", particle_list[i].x, particle_list[i].y,particle_list[i].index);
                }

                fclose(gif);
            #endif
        }

        if(number_of_clusters == 1){
            ++cluster_one_counter;
        }

        //Possible Break conditinos to end while loop
        if(steps_taken >= break_condition){
            printf("Stopped after too many steps: %d\n",steps_taken);
            break;
        }

        if((prob_separate == 0) && (number_of_clusters == 1)){
            printf("1 Cluster Reached\n");
            break;
        }

        if(cluster_one_counter > progress){
            printf("Clusters Oscilating over 1 final cluster.\n");
            break;
        }

        steps_taken++;
        connected = False;
        separation = False;
    }

    free(connecting_clusters);
    free(odist);

    printf("Total number of clusters %d\n", number_of_clusters);
    #ifdef GIF
        sprintf(gif_file_name, "Animation/%d-ClusterSize%dParticles%dProb%f.csv", steps_taken, lat_size, num_particles, prob_separate);
        gif = fopen(gif_file_name, "w");

        for(int i = 0; i < num_particles; i++){
            fprintf(gif, "%lf,%lf,%d\n", particle_list[i].x, particle_list[i].y,particle_list[i].index);
        }

        fclose(gif);
    #endif
    //Delete Partial Results File
    remove(middle_file_name);

}

void saveResults(){
       // Create result directory if directory does not exist
    struct stat r = {0};
    if (stat("/Results", &r) == -1){
        #ifdef __unix__
            mkdir("Results", 0777);
        #else
            mkdir("Results");
        #endif
    }

    char resultrg_file_name[80];
    FILE *frg2;

    #ifdef MAX_COORDINATION
        sprintf(resultrg_file_name, "Results/Rg2EdgeClusterSize%dParticles%dProb%f.txt", lat_size, num_particles, prob_separate);
    #else
        sprintf(resultrg_file_name, "Results/Rg2ClusterSize%dParticles%dProb%f.txt", lat_size, num_particles, prob_separate);
    #endif

    #ifdef CLUSTER_NUMBER
        writeNumberFile();
    #endif

    frg2 = fopen(resultrg_file_name, "w");
    for (int i = 0; i < number_of_clusters; i++)
        fprintf(frg2, "%lf\n", cluster_list[i].rg2);
    fclose(frg2);

    char result_file_name[80];
    FILE *fr;
    #ifdef MAX_COORDINATION
        sprintf(result_file_name, "Results/EdgeClusterSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #else
        sprintf(result_file_name, "Results/ClusterSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #endif

    fr = fopen(result_file_name, "w");
    for(int i = 0; i < num_particles; i++){
        fprintf(fr, "%lf,%lf,%d,%d\n", particle_list[i].x, particle_list[i].y, particle_list[i].index, particle_list[i].coordination_number);
    }
    fclose(fr);
}

//Function to allocate memory to global pointers
void allocate_memory(){
    firstp = (int *)malloc(sizeof(int) * num_particles);
    nextp = (int *)malloc(sizeof(int) * num_particles);
    lastp = (int *)malloc(sizeof(int) * num_particles);
    particle_list = (Particle *)malloc(sizeof(Particle) * num_particles);
    cluster_list = (Cluster *)malloc(sizeof(Cluster) * num_particles);
    cell_list = (int *)malloc(sizeof(int) * ((lat_size * lat_size) + num_particles));
    mass_list = (int *)malloc(sizeof(int) * num_particles);

    k_list = (int **)malloc(sizeof(int *) * (lat_size * lat_size));
    for (int i = 0; i < (lat_size * lat_size); i++){
        k_list[i] = (int *)malloc(sizeof(int) * k_next);
    }
    list_z1.particle_list = (StackCircular **)malloc(sizeof(StackCircular) * num_particles);
}

//Function to deallocate memory to global pointers
void deallocate_memory(){
    free(firstp);
    free(nextp);
    free(lastp);
    free(particle_list);
    free(cluster_list);
    free(cell_list);
    free(mass_list);
    free(k_list);
    free(list_z1.particle_list);
}

// First function ran in the procedure, creates all necesary structures and places all particles in the system
void initialize(){
    int allocated_progress = 10000;
    cells = lat_size; // Set amount of horizontal cells
    cells2 = lat_size * lat_size; // Total number of cells

    // Set all values of cell)_list to -1
    for(int i = 0; i < (cells2 + num_particles); i++) {
        *(cell_list + i) = -1;
    }

    int *kptr;
    kptr = (int *)malloc(sizeof(int) * k_next);
    // Set all combinations of the cells in k
    for(int k = 0; k < cells2; k++){
        setKlist(k, kptr);
    }
    free(kptr);

    //Initialize Z1 list
    list_z1.total = 0;
    list_z1.start = NULL;

    // Set base values outside the system for all particles in system
    // This is done to avoid a clash when comparing elements of not created list in checkSpot.
    // Also sets all values in the mass_list to 0.
    for(int i = 0; i < num_particles; i++){
        Particle particle;
        particle.x = -2 * cells2;
        particle.y = -2 * cells2;
        particle.k = -cells2;
        particle_list[i] = particle;

        mass_list[i] = 0;
    }
    number_of_clusters = 0;

    double current_x, current_y;
    double initial_rg2 = ring_size / 8;
    for(int i = 0; i < num_particles; i++){
        current_x = ((double) rand() / (double)(RAND_MAX)) * lat_size;
        current_y = ((double) rand() / (double)(RAND_MAX)) * lat_size;
        do{
            current_x = ((double) rand() / (double)(RAND_MAX)) * lat_size;
            current_y = ((double) rand() / (double)(RAND_MAX)) * lat_size;
        }while(checkSpot(current_x, current_y) == False);
        // Places the resulting values in the current strcture of Particle
        particle_list[i].x = current_x;
        particle_list[i].y = current_y;
        particle_list[i].number = i;
        particle_list[i].index = number_of_clusters;
        setK(particle_list + i);
        particle_list[i].coordination_number = 0;
        particle_list[i].neighbor = NULL;

        //Create stack node for each particle
        StackCircular *new = (StackCircular *)malloc(sizeof(StackCircular));
        new->number = i;
        new->next = NULL;
        new->prev = NULL;
        list_z1.particle_list[i] = new;


        // Places the particle (Cluster) on the firstp, nextp, and lastp lists
        firstp[number_of_clusters] = i;
        nextp[i] = -1;
        lastp[number_of_clusters] = number_of_clusters;


        Cluster cluster; // Initialize a Cluster struct
        // Allocate corresponding values of mass, center of mass and radius of gyration squared.
        cluster.mass = 1;
        cluster.cx = current_x;
        cluster.cy = current_y;
        // cluster.rg2 = ((double)ring_size) / ((double)8);
        cluster.rg2 = 0;
        cluster_list[i] = cluster;

        addToCellList(particle_list + i); // Adds particle to the cell_list

        number_of_clusters++;

        if (i % allocated_progress == 0)
            printf("%i Allocated Clusters\n", i);
    }
    mass_list[0] = num_particles; // Set mass list[0] to # of particles
    denominator = num_particles; // Sets denominator: denominator = (1/m_c1) + (1/m_c2) + ... + (1/m_cn) = 1 * n
    setA(); // Sets constant of normalization.

    #ifdef RGINFO
        resetRgFile();
    #endif


    #ifdef CLUSTER_NUMBER
        writeNumberFile();
    #endif
    printf("%i Allocated Clusters\n", num_particles);
    printf("\n");
}

// Function which returns the 9 cells that are neighbors of a given k cell.
void setKlist(int k, int *kptr){
    // The set list for every k will be as follows:
    // [ Cell bellow to the left, Cell bellow, Cell bellow to the right,
    //   Cell to the left, Current cell (k), Cell to the right,
    //   Cell above to the left, Cell above, Cell above to the right ]

    if (k == 0){
        int k_values[k_next] = {cells2 - 1, cells2 - cells, cells2 - cells + 1, cells - 1, 
                                k, k + 1, k + (2 * cells) - 1, k + cells, k + cells + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];

    } // Bottom left corner of lattice

    else if (k == cells2 - 1){
        int k_values[k_next] = {k - cells - 1, k - cells, k - (2 * cells) + 1, 
                                k - 1, k, k - cells + 1, cells - 2, cells - 1, 0};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    }  // Top right corner of lattice

    else if (k == cells2 - cells){
        int k_values[k_next] = {k - 1, k - cells, k - cells + 1, cells2 - 1, 
                                k, k + 1, cells - 1, 0, 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // Top left corner of lattice

    else if (k == cells - 1){
        int k_values[k_next] = {cells2 - 2, cells2 - 1, cells2 - cells, k - 1, k, 
                                k - cells + 1, k + cells - 1, k + cells, k + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // Bottom right corner of lattice

    else if (k % lat_size == 0){
        int k_values[k_next] = {k - 1, k - cells, k - cells + 1, k + cells - 1, k,
                                k + 1, k + (2 * cells) - 1, k + cells, k + cells + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // West Border of lattice

    else if (k % lat_size == lat_size - 1){
        int k_values[k_next] = {k - cells - 1, k - cells, k - (2 * lat_size) + 1, k - 1, 
                                k, k - cells + 1, k + cells - 1, k + cells, k + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // East border of lattice

    else if ((int)k / lat_size == 0){
        int k_values[k_next] = {k + cells2 - cells - 1, k + cells2 - cells, 
                                k + cells2 - cells + 1,  k - 1, k, k + 1, 
                                 k + cells - 1, k + cells, k + cells + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // South border of lattice

    else if ((int)k / lat_size == lat_size - 1){
        int k_values[k_next] = {k - cells - 1, k - cells, k - cells + 1, k - 1, 
                                k, k + 1, k - cells2 + cells - 1, k - cells2 + cells,
                                k - cells2 + cells + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // North border of lattice

    else{
        int k_values[k_next] = {k - cells - 1, k - cells, k - cells + 1, k - 1, 
                                k, k + 1, k + cells - 1, k + cells, k + cells + 1};
        // kptr = k_values;
        for(int i = 0; i < k_next; ++i)
            kptr[i] = k_values[i];
    } // Base list for any non border case in the grid

    for(int i = 0; i < k_next; i++){
        k_list[k][i] = kptr[i];
    }  
}

// Adds a particle's k to the Cell List to have a relative position of all the particles
void addToCellList(Particle *particle){
    int past;
    int current = (*particle).k;

    while (current != -1){
        past = current;
        current = cell_list[current];
    } // If that k space is occupied the particle is linked to the original particles in the list

    cell_list[past] = cells2 + (*particle).number;
}

// Removes one particle from the cell list (to later be placed again after moving)
void resetCellListElement(Particle *particle){
    int current;
    int position = (*particle).k;

    current = cell_list[position];
    while (current != -1){
        if (current - cells2 == (*particle).number){
			cell_list[position] = cell_list[current];
			cell_list[current] = -1;
            break;
        } 

        position = current;
        current = cell_list[position];
    } // Runs until the specified particle is found
}

// Removes one particle from the next particle list (since it separates from the cluster)
void resetParticleInLists(int number){
    int current, i = 0;
    int position = particle_list[number].index;

    current =  firstp[position];

    if (number == firstp[position]){
        firstp[position] = nextp[firstp[position]];
        nextp[current] = -1;
    }
	else{
		while (current != -1){
			if (current == number){
                nextp[position] = nextp[current];
		        nextp[current] = -1;
				break;
			} 
			position = current;
			current = nextp[position];
            i++;
		} // Runs until the specified particle is found
		
		if (number == lastp[particle_list[number].index]){
            if (i == 1){
                lastp[particle_list[number].index] = position;
            }
            else{
                lastp[particle_list[number].index] = position;
            }
		}
	}

	firstp[number_of_clusters] = number;
	lastp[number_of_clusters]  = number;
	
}

void removeParticleZ1Lists(int number){
    if (list_z1.particle_list[number] == list_z1.start){
        list_z1.start = list_z1.start->next;
    }

    popList(list_z1.particle_list[number]);
    list_z1.total--;
}

void addParticleZ1List(int number){
    if (list_z1.start != NULL){
        list_z1.start->prev = list_z1.particle_list[number];
    }

    list_z1.particle_list[number]->next = list_z1.start;
    list_z1.particle_list[number]->prev = NULL;
    list_z1.start = list_z1.particle_list[number];
    list_z1.total++;

}

// When joining clusters the masses are joined to change the value of the denominator correspondingly
void changeDenominator(int mass1, int mass2){
    denominator += (1 / (double)(mass1 + mass2)) - (1 / (double)mass1) - (1 / (double)mass2);
    if (denominator <= 0) printf("BREAK"); // Created to print if the value of the denominator is an "impossible value"
}

// When joining clusters the masses are joined and changed in the mass_list correspondingly
void changeMassList(int mass1, int mass2){
    mass_list[mass1 - 1]--;
    mass_list[mass2 - 1]--;
    mass_list[(mass1 + mass2) - 1]++;
}

// With a given denominator the constant of normalization is calculated
void setA(){
    A = 1 / denominator;
}

int selectClusterMass(){
    int selected_cluster;
    double prob_num, prob_den, prob, z;

    // Do while loop to select cluster to move, if the probability of the cluster is higher than z the cluster moves.
    do{
        selected_cluster = rand() % number_of_clusters; // Random cluster index from remaining clusters

        prob_num = mass_list[cluster_list[selected_cluster].mass - 1] * (A);
        prob_den = cluster_list[selected_cluster].mass;
        prob = prob_num /  prob_den; // Probability

        z = (double)rand() / (double)RAND_MAX; // Random number between 0 and 1.
    }while(z > prob);

    return selected_cluster;
}

// Function which takes a position and checks if this position is within interaction range of another particle
int checkSpot(double x, double y){
    double dx, dy, distance;
    int current_k = ((int) (y / ring_size) * cells) + ((int) x / ring_size);

    if (current_k >= cells2 || current_k < 0) 
        return False; 
    // Checks if the particle is outside of the allowed box (special case for y = 500 (exactly))

    for(int i = 0; i < k_next; i++){
        for(int j = 0; j < num_particles; j++){
            if (particle_list[j].k == k_list[current_k][i]){
                dx = fabs(x - (*(particle_list + j)).x);
                dy = fabs(y - (*(particle_list + j)).y);

                if (dx > (lat_size / 2)){
                    dx -= lat_size;
                }

                if (dy > (lat_size / 2)){
                    dy -= lat_size;
                }

                distance = sqrt((dx * dx) + (dy * dy));
                if (distance <= (double)ring_size){
                    return False;
                }
            }
        }
    }
    return True;
}

// Given a cluster and a direction the cluster is moved in a specified direction, as well as changing cell_list accordingly
void step(int selected_cluster, double dir){
    double dir_x = ring_size * cos(dir);
    double dir_y = ring_size * sin(dir);

    int particle;   

    cluster_list[selected_cluster].cx += dir_x;
    cluster_list[selected_cluster].cy += dir_y;

    cluster_list[selected_cluster].cx = periodicBoundaryConditions(cluster_list[selected_cluster].cx);
    cluster_list[selected_cluster].cy = periodicBoundaryConditions(cluster_list[selected_cluster].cy);

    particle = firstp[selected_cluster];
    while (particle != -1){
        resetCellListElement(particle_list + particle);

        particle_list[particle].x += dir_x;
        particle_list[particle].y += dir_y;

        particle_list[particle].x = periodicBoundaryConditions(particle_list[particle].x);
        particle_list[particle].y = periodicBoundaryConditions(particle_list[particle].y);

        setK(particle_list + particle); // Sets the k value of the moved particle.
        addToCellList(particle_list + particle); // Adds it to the cell list.

        particle = nextp[particle];
    }
}

// Given a cluster, a direction, and a return distance the cluster is moved in a specified direction to remove overlap, 
    // as well as changing cell_list accordingly
void stepBack(int selected_cluster, double *odist, double dir){
    double dir_x = ((*odist) - ring_size) * cos(dir);
    double dir_y = ((*odist) - ring_size) * sin(dir);

    int particle;

    cluster_list[selected_cluster].cx += dir_x;
    cluster_list[selected_cluster].cy += dir_y;

    cluster_list[selected_cluster].cx = periodicBoundaryConditions(cluster_list[selected_cluster].cx);
    cluster_list[selected_cluster].cy = periodicBoundaryConditions(cluster_list[selected_cluster].cy);

    particle = firstp[selected_cluster];
    while (particle != -1){
        resetCellListElement(particle_list + particle);
        particle_list[particle].x += dir_x;
        particle_list[particle].y += dir_y;

        particle_list[particle].x = periodicBoundaryConditions(particle_list[particle].x);
        particle_list[particle].y = periodicBoundaryConditions(particle_list[particle].y);

        setK(particle_list + particle); // Sets the k value of the moved particle.
        addToCellList(particle_list + particle); // Adds it to the cell list.

        particle = nextp[particle];
    }
}

// If a particle is placed outside the periodic boundary conditions the function pushes it back on the opposite side.
double periodicBoundaryConditions(double val){
    if (val > lat_size){
        val -= (double)lat_size;
    }
    else if (val < 0){
        val += (double)lat_size;
    }
        
    return val;
}

// Given a cluster all values from *cluster and *odist are set to base values to then be compared, checks for every particle in the cluster
int checkCluster(int selected_cluster, double dir, int *clusters, double *odist){
    int particle;

    clusters[0] = clusters[1] = clusters[2] = clusters[3] = -1; // Initially sets all values to -1 to avoid issues.
    *odist = 2 * ring_size;

    particle = firstp[selected_cluster];
    while (particle != -1){
        checkParticle(particle, dir, clusters, odist);
        particle = nextp[particle];
    }

    if (*odist != 2 * ring_size){
        return True;
    }
    else{
        return False;
    }
}

// Chekcs if a particle from another cluster is in interaction range from it
void checkParticle(int selected_particle, double dir, int *checked_clusters, double *odist){
    int test_part;
    int i;
    double p1x = particle_list[selected_particle].x, p1y = particle_list[selected_particle].y; // Position of selected particle
    double p2x, p2y, dx, dy, dist;

    double new_odist;

    int kpar;
    for(i = 0; i < k_next; i++){
        kpar = k_list[particle_list[selected_particle].k][i];
        while (cell_list[kpar] != -1){
            test_part = cell_list[kpar] - cells2;

            p2x = (particle_list[test_part]).x;
            p2y = (particle_list[test_part]).y;

            dx = p1x - p2x;
            dy = p1y - p2y;

            if (dx > (lat_size / 2))
                dx -= lat_size;
            if (dx < (-lat_size / 2))
                dx += lat_size;

            if (dy > (lat_size / 2))
                dy -= lat_size;
            if (dy < (-lat_size / 2))
                dy += lat_size;

            dist = sqrt((dx * dx) + (dy * dy));
            if ((dist <= (double)ring_size) && (particle_list[selected_particle].index != particle_list[test_part].index)){
                new_odist = overlapDist(dist, dir, p1x, p1y, p2x, p2y); // Calculates distance

                // If stored distance is less than the previos (more overlap) this particle becomes the connecting particle
                if (new_odist < *odist){
                    *(checked_clusters + 0) = particle_list[selected_particle].number;
                    *(checked_clusters + 1) = particle_list[selected_particle].index;
                    *(checked_clusters + 2) = particle_list[test_part].number;
                    *(checked_clusters + 3) = particle_list[test_part].index;
                    *odist = new_odist;
                }

            }
            kpar = cell_list[kpar];

        }
    }
}


// Calculates the overlap distance of two neigboring particles
float overlapDist(double dist, double dir, double p1x, double p1y, double p2x, double p2y){
    double x0 = periodicBoundaryConditions(p1x - (ring_size * cos(dir)));
    double y0 = periodicBoundaryConditions(p1y - (ring_size * sin(dir)));

    double dx = p2x - x0;
    double dy = p2y - y0;
    if (dx > (lat_size / 2))
        dx -= lat_size;
    if (dx < (-lat_size / 2))
        dx += lat_size;

    if (dy > (lat_size / 2))
        dy -= lat_size;
    if (dy < (-lat_size / 2))
        dy += lat_size;

    double x1 = fabs((dx * cos(dir)) + (dy * sin(dir)));

    double new_odist = x1 - sqrt((ring_size * ring_size) - (dist * dist) + ((x1 - ring_size) * (x1 - ring_size)));

    return new_odist;
}

// Function to join two clusters by index
void joinClusters(int c1, int c2){
    int lc_label, sc_label, lc_mass, sc_mass;
    // lc_label = larger cluster label and, sc_label = smaller cluster label
    if (cluster_list[c1].mass > cluster_list[c2].mass){
        lc_label = c1;
        lc_mass = cluster_list[c1].mass;
        sc_label = c2;
        sc_mass = cluster_list[c2].mass;
    }
    else{
        lc_label = c2;
        lc_mass = cluster_list[c2].mass;
        sc_label = c1;
        sc_mass = cluster_list[c1].mass;
    }
        
    // Manipulation of linked lists
    nextp[lastp[lc_label]] = firstp[sc_label];
    lastp[lc_label] = lastp[sc_label];

    changeDenominator(lc_mass, sc_mass);
    changeMassList(lc_mass, sc_mass);

    cluster_list[lc_label].mass += cluster_list[sc_label].mass;

    int current = firstp[sc_label];
    while (current != -1){
        particle_list[current].index = lc_label;
        current = nextp[current];
    }

    number_of_clusters--;

    double cm_x_list[2] = {cluster_list[lc_label].cx, cluster_list[sc_label].cx};
    double cm_y_list[2] = {cluster_list[lc_label].cy, cluster_list[sc_label].cy};
    int cm_mass_list[2] = {lc_mass, sc_mass};

    double *cx = (double *)malloc(sizeof(double)), *cy = (double *)malloc(sizeof(double));

    centerOfMassTwoPoints(cm_x_list, cm_y_list, cm_mass_list, cx, cy);
    
    cluster_list[lc_label].rg2 = radiusOfGyration(lc_mass, sc_mass, lc_label, sc_label, cx, cy);
    cluster_list[lc_label].cx = *cx;
    cluster_list[lc_label].cy = *cy;


    // Resets label of the grid to the smaller of the two labels to have the last cluster existing be a 0 cluster.
    if (sc_label != number_of_clusters){
        int current = firstp[number_of_clusters];
        while (current != -1){
            particle_list[current].index = sc_label;
            current = nextp[current];
        }
        
        // Moves the characteristics of the largest index cluster to a smaller available index
        firstp[sc_label] = firstp[number_of_clusters];
        lastp[sc_label] = lastp[number_of_clusters];
        cluster_list[sc_label].mass = cluster_list[number_of_clusters].mass;
        cluster_list[sc_label].cx = cluster_list[number_of_clusters].cx;
        cluster_list[sc_label].cy = cluster_list[number_of_clusters].cy;
        cluster_list[sc_label].rg2 = cluster_list[number_of_clusters].rg2;
    }

    setA();

    free(cx);
    free(cy);

    #ifdef RGINFO
        if (cluster_list[lc_label].mass > 2)
            writeRgFile(lc_label);
    #endif

}

int separateCluster(int number){
    int adjacent_particle;
	int lc_index = particle_list[number].index;
	int lc_mass = cluster_list[lc_index].mass;

	int sc_mass = 1;
    adjacent_particle = particle_list[number].neighbor->number;

    if(particle_list[number].coordination_number == 1){
        removeParticleZ1Lists(number);
    } 

    if(particle_list[adjacent_particle].coordination_number == 1){
        removeParticleZ1Lists(adjacent_particle);
    }

    --particle_list[number].coordination_number;
    --particle_list[adjacent_particle].coordination_number;

    if(particle_list[adjacent_particle].coordination_number == 1){
        addParticleZ1List(adjacent_particle);
    }

    resetParticleInLists(number);
	
	double cm_x_list[2] = {cluster_list[lc_index].cx, particle_list[number].x};
    double cm_y_list[2] = {cluster_list[lc_index].cy, particle_list[number].y};
    int cm_mass_list[2] = {lc_mass - 1, sc_mass};

    double *new_cx = (double *)malloc(sizeof(double)), *new_cy = (double *)malloc(sizeof(double));
    undoCenterOfMassTwoPoints(cm_x_list, cm_y_list, cm_mass_list, new_cx, new_cy);

    double new_rg2 = undoRadiusOfGyration(sc_mass, lc_index, number, new_cx, new_cy);

	cluster_list[lc_index].cx = *new_cx;
	cluster_list[lc_index].cy = *new_cy;
	cluster_list[lc_index].rg2 = new_rg2;
    --cluster_list[lc_index].mass;

    #ifdef RGINFO
        if (cluster_list[lc_index].mass > 2)
            writeRgFile(lc_index);
    #endif
	
    particle_list[number].index = number_of_clusters;
    cluster_list[number_of_clusters].mass = 1;
    cluster_list[number_of_clusters].cx = particle_list[number].x;
    cluster_list[number_of_clusters].cy = particle_list[number].y;
    cluster_list[number_of_clusters].rg2 = 0;

	denominator += (-1 / (double)lc_mass) + (1 / (double)(lc_mass - 1)) + 1;
	mass_list[lc_mass - 1]--;
	mass_list[lc_mass - 2]++;
	mass_list[0]++;
    setA();

    ++number_of_clusters;

    particle_list[number].neighbor = pop(particle_list[number].neighbor);
	
	particle_list[adjacent_particle].neighbor = popVal(number, particle_list[adjacent_particle].neighbor);    

    free(new_cx);
    free(new_cy);

    return adjacent_particle;
}


// Function that gives the position of all particles in a cluster utilizing the linked list
void posParticlesCluster(int cluster_label, double *x_list, double *y_list){
    int current = firstp[cluster_label];
    int i = 0;
    while (current != -1){
        x_list[i] = particle_list[current].x;
        y_list[i] = particle_list[current].y;

        current = nextp[current];
        i++;
    } 
}

void centerOfMassTwoPoints(double *cx_list, double *cy_list, int *mass_list, double *cx, double *cy){
    int total_mass = mass_list[0] + mass_list[1];
    double dx = cx_list[1] - cx_list[0], dy = cy_list[1] - cy_list[0];

    if (dx > (lat_size / 2))
        dx -= lat_size;
    if (dx < (-lat_size / 2))
        dx += lat_size;

    if (dy > (lat_size / 2))
        dy -= lat_size;
    if (dy < (-lat_size / 2))
        dy += lat_size;

    *cx = periodicBoundaryConditions(cx_list[0] + (((double)mass_list[1] / (double)total_mass) * dx));
    *cy = periodicBoundaryConditions(cy_list[0] + (((double)mass_list[1] / (double)total_mass) * dy));

}

void undoCenterOfMassTwoPoints(double *cx_list, double *cy_list, int *mass_list, double *cx, double *cy){
    int total_mass = mass_list[0] + mass_list[1];
    double dx = cx_list[1] - cx_list[0], dy = cy_list[1] - cy_list[0];

    if (dx > (lat_size / 2))
        dx -= lat_size;
    if (dx < (-lat_size / 2))
        dx += lat_size;

    if (dy > (lat_size / 2))
        dy -= lat_size;
    if (dy < (-lat_size / 2))
        dy += lat_size;

    *cx = periodicBoundaryConditions(cx_list[0] - (((double)mass_list[1] / (double)total_mass) * dx));
    *cy = periodicBoundaryConditions(cy_list[0] - (((double)mass_list[1] / (double)total_mass) * dy));
}

double radiusOfGyration(int lc_mass, int sc_mass, int lc_label, int sc_label, double *cx_new, double *cy_new){
    double dcx1 = fabs((*cx_new) - cluster_list[lc_label].cx), dcy1 = fabs((*cy_new) - cluster_list[lc_label].cy);
    double dcx2 = fabs((*cx_new) - cluster_list[sc_label].cx), dcy2 = fabs((*cy_new) - cluster_list[sc_label].cy);
    if (dcx1 > (lat_size / 2)){
        dcx1 = lat_size - dcx1;
    }
    if (dcy1 > (lat_size / 2)){
        dcy1 = lat_size - dcy1;
    }

    if (dcx2 > (lat_size / 2)){
        dcx2 = lat_size - dcx2;
    }
    if (dcy2 > (lat_size / 2)){
        dcy2 = lat_size - dcy2;
    }

    double h1 = pow(dcx1, 2) + pow(dcy1, 2);
    double h2 = pow(dcx2, 2) + pow(dcy2, 2);

    double rg2_cm11 = cluster_list[lc_label].rg2 + h1;
    double rg2_cm12 = cluster_list[sc_label].rg2 + h2;

    double rg2_f = ((lc_mass * rg2_cm11) + (sc_mass * rg2_cm12)) / (lc_mass + sc_mass);
        
    return rg2_f;
}

double undoRadiusOfGyration(int sc_mass, int lc_label, int removed_particle, double *cx_new, double *cy_new){
    double dcx1 = fabs(cluster_list[lc_label].cx - (*cx_new)), dcy1 = fabs(cluster_list[lc_label].cy - (*cy_new));
    double dcx2 = fabs(cluster_list[lc_label].cx - particle_list[removed_particle].x), dcy2 = fabs(cluster_list[lc_label].cy - particle_list[removed_particle].y);
    if (dcx1 > (lat_size / 2)){
        dcx1 = lat_size - dcx1;
    }
    if (dcy1 > (lat_size / 2)){
        dcy1 = lat_size - dcy1;
    }

    if (dcx2 > (lat_size / 2)){
        dcx2 = lat_size - dcx2;
    }
    if (dcy2 > (lat_size / 2)){
        dcy2 = lat_size - dcy2;
    }

    double h1 = pow(dcx1, 2) + pow(dcy1, 2);
    double h2 = pow(dcx2, 2) + pow(dcy2, 2);

    double rg2_1 = ((double)((cluster_list[lc_label].rg2 * cluster_list[lc_label].mass) - (sc_mass * h2)) / (double)(cluster_list[lc_label].mass - sc_mass)) - h1;
        
    return rg2_1;
}

//Opens RgMassTime file with w so that the file is effectively erased
void resetRgFile(){
    FILE *fp;

    char file_name[60];
    #ifdef MAX_COORDINATION
        sprintf(file_name, "Results/RgMassTimeEdgeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #else
        sprintf(file_name, "Results/RgMassTimeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #endif
    

    fp = fopen(file_name, "w");
    fclose(fp);
}

// Writes a new cluster into the RgMassTime file, this is done every time a cluster joins itself.
void writeRgFile(int cluster_index){
    FILE *fp;

    char file_name[60];
    
    #ifdef MAX_COORDINATION
        sprintf(file_name, "Results/RgMassTimeEdgeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #else
        sprintf(file_name, "Results/RgMassTimeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #endif

    fp = fopen(file_name, "a");
    fprintf(fp, "%d,%lf\n", cluster_list[cluster_index].mass, cluster_list[cluster_index].rg2);

    fclose(fp);

}

void resetNumberFile(){
    FILE *fp;

    char file_name[60];
    #ifdef MAX_COORDINATION
        sprintf(file_name, "Results/NumberTimeEdgeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #else
        sprintf(file_name, "Results/NumberTimeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #endif
    

    fp = fopen(file_name, "w");
    fclose(fp);
}

// Writes a new cluster into the RgMassTime file, this is done every time a cluster joins itself.
void writeNumberFile(){
    FILE *fp;

    char file_name[60];
    
    #ifdef MAX_COORDINATION
        sprintf(file_name, "Results/NumberTimeEdgeSize%dParticles%dProb%f.csv", lat_size, num_particles,prob_separate);
    #else
        sprintf(file_name, "Results/NumberTimeSize%dParticles%dProb%f.csv", lat_size, num_particles, prob_separate);
    #endif

    fp = fopen(file_name, "a");
    fprintf(fp, "%d,", steps_taken);
    for (int i = 0; i < num_particles; ++i){
        fprintf(fp, "%d", mass_list[i]);
        if(i != num_particles - 1)
            fprintf(fp, ", ");
        else
            fprintf(fp, "\n");
    }
    

    fclose(fp);

}

