#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>

#define Pi 3.14159265358979
#define True 1
#define False 0

// Code entries
int num_particles;
int lat_size;
int ring_size;

double *x_list, *y_list;
int *k_list;

int main(int argc, char *argv[]){

	char *file_name = (char *)malloc(sizeof(char) * 100); // Allocate memory to filename

	double time_spent = 0.0;
    clock_t begin = clock(); // Start timer

    if (argc == 4){
        lat_size = atoi(argv[1]);
        num_particles = atoi(argv[2]);
        strcpy(file_name, argv[3]);
    }
    else{
    	printf("Incorrect number of entries.\n");
    }

   	ring_size = 1;

    x_list = (double *)malloc(sizeof(double) * num_particles);
    y_list = (double *)malloc(sizeof(double) * num_particles);
    k_list = (int*)malloc(sizeof(int) * num_particles);

    FILE *cluster;
    cluster = fopen(file_name, "r");
    char line[256];
    free(file_name); // Remove file_name space (it is not usefull once the lists are created.)

    int line_number = 0;
    double x, y;
    while(fgets(line, sizeof(line), cluster)){
    	sscanf(line, "%lf,%lf", &x, &y);
    	
    	k_list[line_number] = ((int) (y / ring_size) * lat_size) + ((int) (x / ring_size)); // Sets k value of a particle struct
        line_number++;
    }
    fclose(cluster);


    int percolation = False, x_percolation = True, y_percolation = True;
    int break_condition_x, break_condition_y;
    int current_kx_left = False, current_ky_bottom = False;
    int i = 0;
    do{
    	if (x_percolation){
    		break_condition_x = False; // Ends one instance of wile and for loop to move on to the next column
            for (int j = 0; j < lat_size; j++){// This exterior for will check in every column space for a row
            	current_kx_left = (j * lat_size) + i;

                for (int k = 0; k < num_particles; k++){// Checks every particle for the given column
                	if (k_list[k] == current_kx_left){
                    	break_condition_x = True; // Breaks loop
                        break;
                    }
                } 
                    
                        
                if (break_condition_x) 
                	break;

            }
                
            if (!break_condition_x) 
            	x_percolation = False;


    	}

        if (y_percolation){
        	break_condition_y = False; // Ends one instance of wile and for loop to move on to the next row
            for (int j = 0; j < lat_size; j++){ // This exterior for will check in every row space for a column
                current_ky_bottom = (i * lat_size) + j;
                
                for (int k = 0; k < num_particles; k++){// Checks every particle for the given column
                	if (k_list[k] == current_ky_bottom){
                    	break_condition_y = True; // Breaks loop
                        break;
                    }
                }

                if (break_condition_y) 
                	break;
             
            }

            if (!break_condition_y) 
            	y_percolation = False;

        }
            
        i++;

        if (i > lat_size){
        	percolation = True;
            break;
        }

    }while((x_percolation == True) || (y_percolation == True));
        
    free(k_list);

    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    printf("Time elpased is %f seconds\n", time_spent);

    return percolation;
}