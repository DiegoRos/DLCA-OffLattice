// ============================================================
// Code: DLCA 2D Fractal Dimension
// File: FracDimDLCA.c
// Author: Diego Rosenberg
/*
Description: This code calculates the fractal dimension of a 2d
object utilizing two techniques. This code returns a double 
representing the calculated fractal dimension.
*/
// ============================================================

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

double *x_list, *y_list;

int divis, break_condition, division_count; 
double *box_counts, *avg_counts, *box_size;


void allocate_memory();
void deallocate_memory();
double fractalDimension(double *x, double *y);
void fractalDimensionCounts(double *x, double *y);
void linearReg(int n, double *x, double *y, double *m, double *b);

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
   	
   	divis = 4 + (4 * (num_particles <= 8000));
	break_condition = 2 + (2 * (lat_size > 100));

	if (break_condition == 4){
		division_count = (int) (log(lat_size / divis) / log(2));
	}
	else{
		division_count = (int) (log(lat_size / (divis / break_condition)) / log(2));
	}


	allocate_memory();
 
    FILE *cluster;
    cluster = fopen(file_name, "r");
    char line[256];
    free(file_name); // Remove file_name space (it is not usefull once the lists are created.)

    int line_number = 0;
    while(fgets(line, 256, cluster)){
    	sscanf(line, "%lf,%lf", (x_list + line_number), (y_list + line_number));
        line_number++;
    }
    fclose(cluster);

    double frac_dim;
	frac_dim = fractalDimension(x_list, y_list);

	struct stat r = {0};
    if (stat("/Results", &r) == -1){
        mkdir("Results", 0777);
    }
    char result_file_name[80];
    FILE *f;
    sprintf(result_file_name, "Results/FracDimCountsSize%dParticles%d.csv", lat_size, num_particles);
    f = fopen(result_file_name, "w");
    for(int i = 0; i < division_count; i++){
        fprintf(f, "%lf,%lf,%lf\n", box_size[i], box_counts[i], avg_counts[i]);
    }
    fclose(f);

   	char fractal_dimension_file_name[80]; 
    sprintf(fractal_dimension_file_name, "Results/FracDimSize%dParticles%d.txt", lat_size, num_particles);
    f = fopen(fractal_dimension_file_name, "w");
    fprintf(f, "%lf\n", frac_dim);
    fclose(f);

	deallocate_memory();

    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    printf("Time elpased is %f seconds\n", time_spent);

	return 0;
}

void allocate_memory(){
	// Allocate Memory to lists:
    x_list = (double *)malloc(sizeof(double) * num_particles);
    y_list = (double *)malloc(sizeof(double) * num_particles);


	// Lists to calculate fractal dimension.
	box_counts = (double *)malloc(sizeof(double) * division_count);
	avg_counts = (double *)malloc(sizeof(double) * division_count);
	box_size = (double *)malloc(sizeof(double) * division_count);
}

void deallocate_memory(){
	free(x_list);
    free(y_list);

    free(box_counts);
	free(avg_counts);
	free(box_size);
}

// Calculates the fractal dimension.
double fractalDimension(double *x, double *y){
	double fractal_dimension;

	fractalDimensionCounts(x, y);

	double *log_box_size = (double *)malloc(sizeof(double) * division_count);
	double *log_box_counts = (double *)malloc(sizeof(double) * division_count);
	double *log_avg_counts = (double *)malloc(sizeof(double) * division_count);

	// Apply logarithm to all count arrays and box_size array.
	for(int i = 0; i < division_count; i++){
		log_box_size[i] = log(box_size[i]);
		log_box_counts[i] = log(box_counts[i]);
		log_avg_counts[i] = log(avg_counts[i]);
	}

	double *m1 = (double *)malloc(sizeof(double)), *b1 = (double *)malloc(sizeof(double));
	double *m2 = (double *)malloc(sizeof(double)), *b2 = (double *)malloc(sizeof(double));
	linearReg(division_count, log_box_size, log_avg_counts, m1, b1);
	linearReg(division_count, log_box_size, log_box_counts, m2, b2);

	fractal_dimension = (fabs(*m1) + fabs(*m2)) / 2;

	return fractal_dimension;
}

// Performs the box count and average count methods to position values
void fractalDimensionCounts(double *x, double *y){
	double division;
	double lower_x, upper_x, lower_y, upper_y;
	int i, j, k;

	int div_layer = 0, box_count; // Counter for how many divisions have been made.
	do{
		division = lat_size / (double)divis;
		box_count = 0;

		for (i = 1; i <= divis; i++){
			for (j = 1; j <= divis; j++){
				lower_x = (i - 1) * division;
				upper_x = i * division;
				lower_y = (j - 1) * division;
				upper_y = j * division;

				for(k = 0; k < num_particles; k++){
					if(( (lower_x <= x[k]) && (x[k] <= upper_x) ) && ( (lower_y <= y[k]) && (y[k] <= upper_y) )){
						box_count++;
						break;
					}
				}
			}
		}

		box_counts[div_layer] = box_count;
		avg_counts[div_layer] = (double)num_particles / (double)box_count;
		box_size[div_layer] = division;

		printf("%lf division size with minimum %d\n", division, break_condition);
		divis = divis * 2;
		div_layer++;
	}while(division >= break_condition);
}

// Calculates linear regression for two lists of doubles of size n and returns the slope and y-intercept in the pointers *m, and *b
void linearReg(int n, double *x, double *y, double *m, double *b){

	double x_sum = 0, y_sum = 0, x_squared_sum = 0, xy_sum = 0;

	for (int i = 0; i < n; i++){
		x_sum += x[i];
		y_sum += y[i];
		x_squared_sum += x[i] * x[i];
		xy_sum += x[i] * y[i];
	}

	*m = ((n * xy_sum) - (x_sum * y_sum)) / ((n * x_squared_sum) - (x_sum * x_sum));
	*b = (y_sum - ((*m) * x_sum)) / n;
}
