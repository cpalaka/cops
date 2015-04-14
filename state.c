#include <stdlib.h>
#include "state.h"
#include "binhash.h"
#include <assert.h>

//function to create space for new sim_state_t
sim_state_t* alloc_state(int n) {
	sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
	s->n = n;
	s->part = NULL;
	s->hash = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));
	return s;
}

void free_state(sim_state_t* s) {
	free(s->hash);
	free(s->part);//CHANGE THIS TO REMOVE ENTIRE LINKED LIST OF PARTICLES. RIGHT NOW IT IS FOR THE OLD DYNAMIC ARRAY HOLDING ALL PARTICLES
	free(s);
}

void remove_particle(int n, sim_state_t* s) {
	particle_t* p = s->part;
	particle_t* prev = s->part;
	assert(n <= s->n && n >0);

	while(p != NULL) {
		if(p->pno == n) {
			if(p == prev) {//head is the element to remove
				prev=prev->ll_next;
				s->part = prev;
				p->ll_next=NULL;
				free(p);
				break;
			}
			if(n == 1) {//last element is the element to remove
				prev->ll_next = NULL;
				free(p);
				break;
			}
			prev->ll_next = p->ll_next;
			p->ll_next=NULL;
			free(p);
			break;
		}
		p->pno = p->pno - 1;
		prev = p;
		p=p->ll_next;
	}

}