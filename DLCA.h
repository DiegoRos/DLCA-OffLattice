// Header file including relevant libraries and structures.
 
#include <stdio.h>
#include <sys/stat.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define Pi 3.14159265358979
#define True 1
#define False 0
// #define MAX_COORDINATION 2

// #define GIF 1
// #define RGINFO 1

typedef struct {
    double x, y; // Particle position
    int number, index; // Particle number (I.D.) and index (cluster I.D.)
    int coordination_number; // Number of particles connected to it
    int k; // Cell position of particle
}Particle;

typedef struct {
    int mass; // Mass of cluster
    double cx, cy; // Center of mass of cluster
    double rg2; // Radius of Gyration squared of cluster
}Cluster;


// Stack struct with push and pop functions.
typedef struct Stack{
	int number; // Saves a particle number
	struct Stack *next; // Saves the next value
}Stack;

Stack * push(int val, struct Stack *s){
	Stack *aux = (Stack *)malloc(sizeof(Stack));
	
	aux->number = val;
	aux->next = s;
	s = aux;
	return s;
}

Stack * pop(Stack *s){
	Stack *aux;
	if (s != NULL){
		aux = s;
		s = s->next;
		free(aux);
	}
	return s;
}

int imprimirStack(Stack *s){
	while (s != NULL){
		printf("Num: %i\n", s->number);
		s = s->next;
	}

	return 0;
}

int get(Stack *s, int pos){
    Stack *aux = s;
	int i = 0;
	if (s == NULL){
		printf("Error gets function recieved NULL stack\n");
		exit(1);
	}

    while (i < pos){
		aux = aux->next;
		++i;
	}
    return aux->number;
}