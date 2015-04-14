#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <time.h>

#include "params.h"
#include "state.h"
#include "binhash.h"
#include "interact.h"
//#include "io.h"

sim_state_t* init_particles(sim_param_t* param) {
	sim_state_t* s = alloc_state(param->N);

	int i, noftypes = param->n_t,ind=0;
	int* netype = (int*) calloc(noftypes, sizeof(int));
	memcpy(netype, param->n_et, sizeof(int)*noftypes);

	//parse particleconfig.txt
	int* initialparticles;
	int* numberifFixed;

	FILE* file = fopen("particleconfig.txt", "r");
	int ch,count=0;
	//count number of lines
    do {
   		ch = fgetc(file);
   		if( ch== '\n') count++;
   	} while( ch != EOF );
   	rewind(file);//need to go back to beginning of file after getting count

   	//holds the initial particle types
   	initialparticles = (int*) calloc(count,sizeof(int));
   	numberifFixed = (int*) calloc(count,sizeof(int));
   	//init all to zero
   	for(i=0; i<count;++i) {
   		numberifFixed[i] = 0;//signifies that its random
   	}

	char line[256];
	int lineno = 0;

	while(fgets(line, sizeof(line), file)) {
		int pflag = 0; //true/1 = fixed, false/0 = random

		for(i=0; i<strlen(line);++i) {
			if(line[i] == ':' || line[i] == ' ') {

				continue;
			}
			if(line[i] == 'f') {
				pflag = 1;
				i+=7;//change place to number

				//get number
				if(pflag) {
					int c = 0;
					char sub[10];//holds number
					int length = strlen(line) - i;
					while (c < length) {
      					sub[c] = line[i+c-1];
      					c++;
   					}
   					sub[c] = '\n';
   					numberifFixed[lineno] = atoi(sub);
   					break;
   				}
			}
			if(i==0) {
				initialparticles[lineno] = line[0] - 97;
			}
		}
		lineno++;
	}
	//find number of total fixed particles (particles required to be read from text file)
	int noffixed = 0;
	for(i=0; i<count;++i) {
		noffixed+=numberifFixed[i];
	}
	
	//start linked list of particles with noffixed particles read from file
	FILE* posFile = fopen("fixedparticles.txt","r");
	int _ind=0;
	for(i = 0; i < noffixed; i++)
  	{
   		float x,y,z;
      	fscanf(posFile,"%f",&x);
      	fscanf(posFile,"%f",&y);
      	fscanf(posFile,"%f",&z);
      	
      	particle_t* particle = (particle_t*) calloc(1, sizeof(particle_t));
      	particle->x[0] = x;
      	particle->x[1] = y;
      	particle->x[2] = z;

      	if(numberifFixed[_ind] > 0) {
			particle->type = _ind;
			numberifFixed[_ind]--;
			netype[_ind]--;
		} else {
			if(_ind == param->n_t -1 && numberifFixed[_ind] == 0) {
				//special case for last particle iterated over
				particle->type = _ind;
			} else {
				particle->type = ++_ind;
				numberifFixed[_ind]--;
				netype[_ind]--;
			}
		}

		particle->ll_next = s->part;
		s->part = particle;
		particle->pno = i+1;
  	}


  	//adjust netype array for already added particles
  	int total_rand_part = param->N - noffixed;


	//random
	for(i=0; i < total_rand_part; ++i) {
		particle_t* newpart = (particle_t*) calloc(1, sizeof(particle_t));

	//randomly give particle a 3d position
		newpart->x[0] =((float)rand()/(float)(RAND_MAX))*param->dim_size;
		newpart->x[1] =((float)rand()/(float)(RAND_MAX))*param->dim_size;
		newpart->x[2] =((float)rand()/(float)(RAND_MAX))*param->dim_size;

		//assign the type of the particle
		if(netype[ind] > 0) {
			newpart->type = ind;
			netype[ind]--;
		} else {
			if(ind == param->n_t -1 && netype[ind] == 0) {
				//special case for last particle iterated over
				newpart->type = ind;
			} else {
				ind++;
				if(netype[ind] != 0) {
					newpart->type = ind;
					netype[ind]--;
				} else {
					ind++;
					newpart->type=ind;
					netype[ind]--;
				}
			}
		}

		//printf("[%f,%f,%f]\n",newpart->x[0],newpart->x[1],newpart->x[2]);
		//add particle to head of particle linked list
		newpart->ll_next = s->part;
		s->part = newpart;
		newpart->pno = noffixed + i+1;
	}

	current_no_of_particles = param->N;
	//hash the particles vv
	hash_particles(s);
	free(netype);
	free(initialparticles);
	free(numberifFixed);
	return s;
}

void write_to_file(FILE* fp, sim_state_t* s) {
	particle_t* p = s->part;
	int i;
	while(p != NULL) {
		fprintf(fp, "p%d: %f %f %f type:%d\n", p->pno, p->x[0], p->x[1], p->x[2], p->type);
		p=p->ll_next;
	}
}

void printpart(sim_state_t* s) {
	particle_t* p = s->part;
	while(p != NULL) {
		printf("p%d: %f %f %f type:%d\n", p->pno, p->x[0], p->x[1], p->x[2], p->type);
		p=p->ll_next;
	}
}
void printn_et(sim_param_t* p) {
	int i;
	for(i=0;i<p->n_t;++i) {
		printf("t%d:%d\n",i,p->n_et[i]);
	}
}

//checks inputs of all reactions and makes sure number of particles is not zero (so that simulation wont fail)
//if they are zero, then terminate simulation
int checkreactions(sim_param_t* p) {
	reaction* r = p->reactions;
	int n = p->M_uni + p->M_bi;

	int i;
	for(i=0; i<n; ++i) {
		if(r[i].type == 1) {//uni
			if(p->n_et[r[i].input[0].ptype] == 0) {
				return -1;
			}
		} else {//bi
			if(p->n_et[r[i].input[0].ptype] == 0) {
				return -1;	
			}
			if(p->n_et[r[i].input[1].ptype] == 0) {
				return -1;
			}
		}
		
	}
	return 1;
}

int main(int argc, char** argv) {

	srand(time(NULL));

	int no_of_timeunits = 6;


	sim_param_t params;
	default_params(&params);
	get_reactions(&params); //get the reactions from reactions.txt

	FILE* fp = fopen(params.fname, "w");
	sim_state_t* state = init_particles(&params);

	int i;
	
	//time loop
printpart(state);
	int tu;
	for(tu = 0; tu < no_of_timeunits; ++tu) {
		int x = checkreactions(&params);
		if(x == -1) {
			printf("Not enough particles for reaction input. Terminating simulation.\n\n");
			break;
		}
		//find reaction with greatest propensity
		reaction* r = compute_propensities(state, &params);
		fire_reaction_event(r,state,&params);

		fire_diffusion_event(state,&params);

		printf("--------------------\n");
		printf("%d ITERATION OVER\n",tu+1);
		printf("--------------------\n");

		free_reaction_list(&params);

	}
	write_to_file(fp, state);

	free_state(state);
	return 0;
}
