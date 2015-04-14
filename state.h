#ifndef STATE_H
#define STATE_H


typedef struct particle_t {
	float x[3]; //particle position
	int type; //particle type
	struct particle_t* next;//used to link particles in same hash bucket
	struct particle_t* ll_next;//used to link all particles in a linked list
	int pno;//particles number or id
} particle_t;

typedef struct sim_state_t {
	int n;            //number of particles
	particle_t* part; //head of particles linked list
	particle_t** hash;//hash table of particles
} sim_state_t;

int current_no_of_particles;
sim_state_t* alloc_state(int n);
void free_state(sim_state_t* s);
void remove_particle(int n, sim_state_t* s);//remove particle with pno=n from the linked list

#endif