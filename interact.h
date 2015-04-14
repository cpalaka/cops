#ifndef INTERACT_H
#define INTERACT_H
#include "params.h"
#include "state.h"
#include "binhash.h"

#define PI 3.14159265

particle_t*** chosen_bi;
particle_t*** chosen_uni;
int* counts_bi;//to keep track of size of dynamically allocated array of each of chosen_bi and chosen_uni
int* counts_uni;
int rtf_count;//size of reaction to fire array
particle_t** reaction_to_fire;

void init_reaction_list(sim_param_t* p);
void free_reaction_list(sim_param_t* p);

float calculate_propensity(int rno, sim_state_t* s, sim_param_t* p);
reaction* compute_propensities(sim_state_t* s, sim_param_t* p);
void fire_reaction_event(reaction* r, sim_state_t* s, sim_param_t* param);
void fire_diffusion_evet(sim_state_t* t, sim_param_t* p);

#endif