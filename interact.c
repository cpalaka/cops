#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <time.h>
#include "interact.h"


//calculate distance between 2 points
static float dist(float* a, float* b) {
	float dx = a[0] - b[0];
	float dy = a[1] - b[1];
	float dz = a[2] - b[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

void init_reaction_list(sim_param_t* p) {

	chosen_uni = (particle_t***) calloc(p->M_uni, sizeof(particle_t**));
	chosen_bi = (particle_t***) calloc(p->M_bi, sizeof(particle_t**));
	counts_uni = (int*) calloc(p->M_uni, sizeof(int));
	counts_bi = (int*) calloc(p->M_bi, sizeof(int));

}

void free_reaction_list(sim_param_t* p) {

	//free chosen uni
	int i;
	for(i=0; i<p->M_uni; ++i) {
		free(chosen_uni[i]);
	}
	free(chosen_uni);

	//free chosen bi
	for(i=0;i<p->M_bi; ++i) {
		free(chosen_bi[i]);
	}
	free(chosen_bi);

	free(counts_uni);
	free(counts_bi);

	for(i=0 ; i<rtf_count; ++i) {
		reaction_to_fire[i] = NULL;
	}
}

float calculate_propensity(int rno, sim_state_t* s, sim_param_t* p) {

	reaction* r = p->reactions;//points to current reaction
	reaction* q = p->reactions;//points to beginning of list of reactions
	r += rno;

	int i;
	int bi_index=0, uni_index=0;//holds uni reaction index and bi reaction index

	particle_t* const* hash = s->hash;
	int cap_rad = p->cap_r;

	float propensity=0;

	if(r->type == 1) {				//xa -> yb

		for(i=0;i<p->M_uni + p->M_bi;++i) {
			if(q == r) break;
			if(r->type == 1) uni_index++;
			q++;
		}

		int a = r->input[0].ptype;
		int x = r->input[0].n;
		x--;//decrement one because we want (x-1) other particles of type a (we exclude the particle we are searching with)

		int p_count = 0;//pcount will keep count of numbers of valid particles
		if(x == 0) {//if x=0, then all particles are considered
			propensity = p->n_et[a] * p->rateconstants[uni_index];

			//choose random particle in set of all particles of type a
			int rndom = (rand()%p->n_et[a])/2;//divide by 2 so it will be faster to find the particle. 

			particle_t** valid_part = (particle_t**) calloc(1, sizeof(particle_t*));
			int acount = 0;
			particle_t* part = s->part;
			while(part != NULL) {
				if(part->type == a) {
					if(acount >= rndom) { //if we find a particle of type a, whose count in all the a particles is greater than the random value we selected, then choose that particle and stop searching
						valid_part[0] = part;
						break;
					}
					acount++;
				}
				part=part->ll_next;
			}
			chosen_uni[uni_index] = valid_part;
			counts_uni[uni_index] = 1;

		} else {//if not, then we need to find how many particles of a are within capture radius of (x-1) other particles of a
			int pcount = 0;
			particle_t* pi = s->part;
			while(pi != NULL) {
				if(pi->type == a) {
					unsigned buckets[MAX_NBR_BINS];
					unsigned noBuckets = particle_neighborhood(buckets, pi);

					particle_t** valid_part = (particle_t**) calloc(x+2, sizeof(particle_t*));
					int valid_count = 0;

					int acount = 0, acheck = 0;
					int j;
					for(j=0;j<noBuckets;++j) {
						particle_t* pj = hash[buckets[j]];
						while(pj != NULL && pi != pj) {
							if(dist(pi->x, pj->x) < cap_rad) {
								if(pj->type == a && acheck==0) {
									if(acount >= x) acheck = 1;
									valid_part[valid_count] = pj;
									valid_count++;
									acount++;
								}
							}
							pj = pj->next;
						}
					}
					if(acheck) {
						//havent randomized selection. for now, the chosen reaction to fire is the last valid one
						valid_part[valid_count] = pi; //add particle used to search (pi) to the very end of list
						chosen_uni[uni_index] = valid_part;
						counts_uni[uni_index] = valid_count;
						pcount++;
					} else {
						free(valid_part);
					}
				}
				pi=pi->ll_next;
			}
			//printf("r(%d) pcount=%d\n",rno,pcount);
			propensity = pcount * p->rateconstants[uni_index];

		}
	} else {							//xa + yb -> zc

		for(i=0;i<p->M_uni + p->M_bi;++i) {
			if(q == r) break;
			if(r->type == 2) bi_index++;
			q++;
		}

		int a = r->input[0].ptype;
		int x = r->input[0].n;
		int b = r->input[1].ptype;
		int y = r->input[1].n;

		int c;//lowest number of a or b (so we can iterate over a lower number of particles)
		int d;//and the other

		int p_count = 0;//count of valid particle pairs

		if(p->n_et[a] < p->n_et[b]) {
			c=a;
			d=b;
			x--;
		}
		else {
			c=b;
			d=a;
			y--;
		}

		int randomselect = rand()%(p->n_et[c]/2);//randomly get a number after which the first valid particle list will be chosen for firing
		int typecnt = 0;//keeps count of the number of molecules of the type which we are iterating over
		int isEventChosen = 0;

		//loop over all particles
		particle_t* pi = s->part;
		while(pi != NULL) {
			if(pi->type == c) {//check only particles of the specified type (in the reaction)
				unsigned buckets[MAX_NBR_BINS];
				unsigned noBuckets = particle_neighborhood(buckets, pi);

				particle_t** valid_part = (particle_t**) calloc(x+y+2, sizeof(particle_t*));//container to hold all valid particles pointers
				int validcount = 0;//count of all  neighbour molecules within capture radius

				int ccount = 0;
				int dcount = 0;//acount and bcount will keep count of valid particles in capture radius(of type a and b)

				int ccheck, dcheck;
				if(c==a) {
					ccheck = (x==0)? 1:0;// acheck and bcheck flags are to stop checking for more particles
					dcheck = (y==0)? 1:0;
				} else {
					ccheck = (y==0)? 1:0;// acheck and bcheck flags are to stop checking for more particles
					dcheck = (x==0)? 1:0;
				}

				int j;
				for(j=0; j<noBuckets; ++j) {//iterate over neighbor buckets
					particle_t* pj = hash[buckets[j]];
					while(pj != NULL && pi != pj) {    //iterate over particles in bucket
						if(dist(pi->x, pj->x) < cap_rad) {  //consider only the ones which are within capture radius

							if(pj->type == c && ccheck==0) {//find x more particles of same type
								if(c==a) {
									if(ccount >= x) ccheck = 1;
								} else {
									if(ccount >= y) ccheck = 1;
								}
								valid_part[validcount] = pj;//add particle to valid list
								validcount++;
								ccount++;
							}

							if(pj->type == d && dcheck == 0) {//make sure the particle is of required type
								if(c==a) {
									if(dcount >= y) dcheck = 1;
								} else {
									if(dcount >= x) dcheck = 1;
								}
								valid_part[validcount] = pj;//add particle to valid list
								validcount++;
								dcount++;
							}
						}
						pj = pj->next;
					}
				}
				//validcount--;//to counteract the increment of validcount after the last valid particle found
				if(ccheck && dcheck) { //if capture radius of pi contains the necessary particles for the reaction, then increment p_count
					if(/*typecnt > randomselect &&*/ !isEventChosen) {//CHECK FOR ERROR HERE: there might not be any particles after the randomselect point
						valid_part[validcount] = pi; // add the particle we used to search with at the last position of array
						chosen_bi[bi_index] = valid_part;
						counts_bi[bi_index] = validcount;
						isEventChosen = 1;

					}
					p_count++;
					//continue;
				} else {
					free(valid_part);
				}
				typecnt++;
			}
			pi=pi->ll_next;
		}
		//printf("p_count:%d\n", p_count);
		propensity = p_count * p->rateconstants[p->M_uni + bi_index]; // propensity = number of valid molecule pairs * rate constant
	}
	return propensity;
}

reaction* compute_propensities(sim_state_t* s, sim_param_t* p) {

	init_reaction_list(p);
	reaction* rs = p->reactions;
	reaction* q = p->reactions;
	int i, numreactions = p->M_uni + p->M_bi;
	//create new array to compare propensities of all reactions
	float* a = (float*) calloc(numreactions, sizeof(float));
	float max =0;
	int result;


	//iterate through all reactions and find maximum
	for(i = 0; i<numreactions; ++i) {

		a[i] = calculate_propensity(i, s, p);

		if(a[i] > max) {
			max = a[i];
			result = i;
		}
	}

	/* tests
	for(i=0; i<numreactions; ++i) {
		printf("propensity r(%d):%f\n\n", i, a[i]);
	}
	//printf("result:%d\n",result);
	*/

	rs += result;

	//set reaction_to_fire for the firing function
	int index = 0;
	if(rs->type == 1) {
		for(i=0;i<numreactions;++i) {
			if(q == rs) break;
			if(q->type==1)index++;
			q++;
		}

		reaction_to_fire = chosen_uni[index];
		rtf_count = counts_uni[index];
	} else {
		for(i=0;i<numreactions;++i) {
			if(q == rs) break;
			if(q->type==2) index++;
			q++;
		}

		reaction_to_fire = chosen_bi[index];
		rtf_count = counts_bi[index];
	}

	free(a);
	return rs;
}

void fire_reaction_event(reaction* r, sim_state_t* state, sim_param_t* param) {
	if(r->type == 1){//unimolecular (xa -> yb)
		int rand_p = rand()%rtf_count;//choose a random particle of the x particles considered to generate an output b within its capture radius
		float rand_p_x[3];
		//remove all particles from reaction_to_fire from the particle linked list
		int i;

		for(i=0;i<rtf_count;++i) {
			if(i==rand_p) {//rand_p_x is the position of the particle which we use to randomly put the output particle around(within capture radius)
				rand_p_x[0] = reaction_to_fire[i]->x[0];
				rand_p_x[1] = reaction_to_fire[i]->x[1];
				rand_p_x[2] = reaction_to_fire[i]->x[2];
			}

			param->n_et[reaction_to_fire[i]->type]--;
			remove_particle(reaction_to_fire[i]->pno, state);

		}

		//update current_no_of_particles
		current_no_of_particles-=rtf_count;

		//add output particle/s to the linked list
		int no_of_outputs = r->output.n;
		for(i=0; i<no_of_outputs; ++i) {
			particle_t* p = (particle_t*) calloc(1, sizeof(particle_t));

			float rand_x = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
			float rand_y = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
			float rand_z = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
			//add the new particle in the capture radius around the previously chosen input particle
			IF_FAIL:if( (rand_p_x[0] + rand_x <= param->dim_size) && (rand_p_x[1] + rand_y <= param->dim_size) && (rand_p_x[2] + rand_z <= param->dim_size) &&
				        (rand_p_x[0] + rand_x >= 0) && (rand_p_x[1] + rand_y >= 0) && (rand_p_x[2] + rand_z >= 0)) {
				p->x[0] = rand_p_x[0] + rand_x;
				p->x[1] = rand_p_x[1] + rand_y;
				p->x[2] = rand_p_x[2] + rand_z;
			} else {
				rand_x = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				rand_y = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				rand_z = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				goto IF_FAIL;
			}


			p->type = r->output.ptype;

			p->ll_next = state->part;
			state->part = p;

			p->pno = ++current_no_of_particles;
			param->n_et[p->type]++;
			//check that none of the particles has left dimension area
			assert(p->x[0] <= param->dim_size && p->x[1] <= param->dim_size && p->x[2] <= param->dim_size);
			assert(p->x[0] >= 0 && p->x[1] >= 0 && p->x[2] >= 0);
		}

	} else { //bimolecular (xa+yb->zc)
		int rand_p = rand()%rtf_count;
		float rand_p_x[3];

		//remove input particles
		int i;
		for(i=0; i<rtf_count; ++i) {
			if(i == rand_p) {
				rand_p_x[0] = reaction_to_fire[i]->x[0];
				rand_p_x[1] = reaction_to_fire[i]->x[1];
				rand_p_x[2] = reaction_to_fire[i]->x[2];
			}
			param->n_et[reaction_to_fire[i]->type]--;
			remove_particle(reaction_to_fire[i]->pno,state);
		}

		//update current_no_of_particles
		current_no_of_particles-=rtf_count;

		//add output particles
		int no_of_outputs = r->output.n;
		for(i=0; i<no_of_outputs;++i) {
			particle_t* p = (particle_t*) calloc(1, sizeof(particle_t));
			
			float rand_x = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
			float rand_y = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
			float rand_z = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;

			//add particle position
			IF_FAIL2:if( (rand_p_x[0] + rand_x <= param->dim_size) && (rand_p_x[1] + rand_y <= param->dim_size) && (rand_p_x[2] + rand_z <= param->dim_size) && 
						 (rand_p_x[0] + rand_x >= 0) && (rand_p_x[1] + rand_y >= 0) && (rand_p_x[2] + rand_z >= 0)) {
				p->x[0] = rand_p_x[0] + rand_x;
				p->x[1] = rand_p_x[1] + rand_y;
				p->x[2] = rand_p_x[2] + rand_z;
			} else {
				rand_x = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				rand_y = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				rand_z = (-1 + ((float)rand()/(float)(RAND_MAX)) + ((float)rand()/(float)(RAND_MAX)))*param->cap_r;
				goto IF_FAIL2;
			}

			p->type = r->output.ptype;

			p->ll_next = state->part;
			state->part = p;

			p->pno = ++current_no_of_particles;
			param->n_et[p->type]++;
			assert(p->x[0] <= param->dim_size && p->x[1] <= param->dim_size && p->x[2] <= param->dim_size);
			assert(p->x[0] >= 0 && p->x[1] >= 0 && p->x[2] >= 0);
		}

	}

	//rehash particles
	hash_particles(state);
}

void fire_diffusion_event(sim_state_t* s, sim_param_t* p) {
	int i;
	float temp = 0;
	int type =0;

	for(i=0;i<p->n_t;++i) {
		float n = p->dconst[i] * p->n_et[i];//propensity of diffusion event( diffusion constant * number of specific molecule type)
		if(temp <= n) {
			type = i;
		}
		temp = n;
	}


	//find random particle of selected type
	particle_t* chosen_part = s->part;

	int rand_part = rand()%p->n_et[type];
	int count=1;
	while(chosen_part != NULL) {
		if(chosen_part->type == type) {
			if(count >= rand_part) break;
			count++;
		}
		chosen_part=chosen_part->ll_next;
	}

	//update the position by randomly moving it in a direction by a constant distance (diff_dist_h)

	float theta = (((float)rand()/(float)(RAND_MAX))*360);//r,theta and phi are used as sperical coordinates for findng the point
	float phi = (((float)rand()/(float)(RAND_MAX))*360);

	theta*=(PI/180);
	phi*=(PI/180);

	float r = p->diff_dist_h[type];//distance to be moved



	float old_x = chosen_part->x[0];
	float old_y = chosen_part->x[1];
	float old_z = chosen_part->x[2];

	//float old[3] = {old_x,old_y,old_z};
	//printf("Particle is at location: (%d,%d,%d).\n", chosen_part->x[0], chosen_part->x[1], chosen_part->x[2]);

	chosen_part->x[0] = (old_x + r*sin(phi)*cos(theta));
	chosen_part->x[1] = (old_y + r*sin(phi)*sin(theta));
	chosen_part->x[2] = (old_z + r*cos(phi));
	//printf("Particle diffused to location: (%d,%d,%d).\n", chosen_part->x[0], chosen_part->x[1], chosen_part->x[2]);

	//assert that the new particle position is inside the limits of max dimension size
	if(!(chosen_part->x[0] <= p->dim_size && chosen_part->x[1] <= p->dim_size && chosen_part->x[2] <= p->dim_size &&
		 chosen_part->x[0] >= 0 && chosen_part->x[1] >= 0 && chosen_part->x[2] >= 0)) {
		chosen_part->x[0] = old_x;
		chosen_part->x[1] = old_y;
		chosen_part->x[2] = old_z;
		fire_diffusion_event(s,p);
	}

	assert(chosen_part->x[0] <= p->dim_size && chosen_part->x[1] <= p->dim_size && chosen_part->x[2] <= p->dim_size);
	assert(chosen_part->x[0] >= 0 && chosen_part->x[1] >= 0 && chosen_part->x[2] >= 0);

	//float distance = dist(chosen_part->x,old);
	//printf("Particle diffuse distance : %f\n", distance);
	//rehash particles
	hash_particles(s);
}