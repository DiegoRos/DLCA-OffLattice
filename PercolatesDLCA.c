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

int *index_list, *k_list;

int largest(int arr[], int n);

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

    index_list = (int *)malloc(sizeof(int) * num_particles);
    k_list = (int *)malloc(sizeof(int) * num_particles);

    printf("%s\n", file_name);

    FILE *fp;
    fp = fopen(file_name, "r");
    
    if(fp == NULL){
        printf("File could not be opened\n");
        exit(1);
    }

    int line_number = 0;
    double x, y;
    int coordination_number;
    while(fscanf(fp, "%lf,%lf,%d,%d\n", &x, &y, &index_list[line_number], &coordination_number) > 1){
    	
    	k_list[line_number] = ((int) (y / ring_size) * lat_size) + ((int) (x / ring_size)); // Sets k value of a particle struct
        line_number++;
    }
    free(file_name); // Remove file_name space (it is not usefull once the lists are created.)
    fclose(fp);

    int num_clusters = largest(index_list, num_particles);

    int percolation, x_percolation, y_percolation, break_condition_x, break_condition_y, current_kx_left, current_ky_bottom;
    int i;

    for(int num = 0; num <= num_clusters; ++num){
        percolation = False, x_percolation = True, y_percolation = True;
        break_condition_x, break_condition_y;
        current_kx_left = False, current_ky_bottom = False;
        i = 0;
        do{
        	if (x_percolation){
        		break_condition_x = False; // Ends one instance of wile and for loop to move on to the next column
                for (int j = 0; j < lat_size; j++){// This exterior for will check in every column space for a row
                	current_kx_left = (j * lat_size) + i;

                    for (int k = 0; k < num_particles; k++){// Checks every particle for the given column
                    	if ((k_list[k] == current_kx_left) && (index_list[k] == num)){
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
                    	if ((k_list[k] == current_ky_bottom)  && (index_list[k] == num)){
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

    if((x_percolation == True) || (y_percolation == True))
        break;
    }
        
    free(k_list);

    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    printf("Time elpased is %f seconds\n", time_spent);

    return percolation;
}

int largest(int arr[], int n){
    int i;
    
    // Initialize maximum element
    int max = arr[0];
 
    // Traverse array elements from second and
    // compare every element with current max 
    for (i = 1; i < n; i++){
        if (arr[i] > max){
            max = arr[i];
        }
    }
 
    return max;
}