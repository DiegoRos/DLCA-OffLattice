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
// #define CLUSTER_NUMBER 1

// Stack struct with push and pop functions.
typedef struct Stack{
	int number; // Saves a particle number
	struct Stack *next; // Saves the next value
}Stack;

Stack * push(int val, Stack *s){
	Stack *aux = (Stack *)malloc(sizeof(Stack));

	if(aux == NULL){
		printf("Not enough memory\n");
		exit(0);
	}
	
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

Stack * popVal(int val, Stack *s){
	Stack *aux = s->next, *prev = s;

	if(s->number == val){
		s = pop(s);
	}
	else{
		while(aux != NULL){
			if(aux->number == val){
				prev->next = aux->next;
				free(aux);
				break;
			}
			prev = aux;
			aux = aux->next;
		}
	}

	return s;
}


// Circular Stack struct with push and pop functions.
typedef struct StackCircular{
	int number; // Saves a particle number
	struct StackCircular *next, *prev; // Saves the next and previous value
}StackCircular;

StackCircular * pushCircular(int val, StackCircular *s){
	StackCircular *aux = (StackCircular *)malloc(sizeof(StackCircular));

	if(aux == NULL){
		printf("Not enough memory\n");
		exit(0);
	}
	aux->number = val;

	if(s != NULL){
		aux->next = s;
		aux->prev = NULL;
		s->prev = aux;
	}
	else{
		aux->next = NULL;
		aux->prev = NULL;
	}
	
	s = aux;
	return s;
}

StackCircular * popList(StackCircular *s){
	StackCircular *aux = s;
	if ((aux->prev == NULL) && (aux->next == NULL)){
		return NULL;
	}
	else if (aux->prev == NULL){
		aux->next->prev = s->prev;
	}
	else if (aux->next == NULL){
		aux->prev->next = s->next;
	}
	else{
		aux->prev->next = s->next;
		aux->next->prev = s->prev;
	}

	s->next = NULL;
	s->prev = NULL;
	return s;
}

int get(StackCircular *s, int pos){
    StackCircular *aux = s;
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

typedef struct{
	int total;
	StackCircular *start;
	StackCircular **particle_list;
}Z1;

typedef struct {
    double x, y; // Particle position
    int number, index; // Particle number (I.D.) and index (cluster I.D.)
    int coordination_number; // Number of particles connected to it
    int k; // Cell position of particle
	Stack *neighbor;
}Particle;

typedef struct {
    int mass; // Mass of cluster
    double cx, cy; // Center of mass of cluster
    double rg2; // Radius of Gyration squared of cluster
}Cluster;
