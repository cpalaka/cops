#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "zmorton.h"
#include "binhash.h"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 * 
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 * 
 *@c*/

#define HASH_MASK (HASH_DIM-1)

void getIndices(particle_t* const p, unsigned* ix, unsigned* iy, unsigned* iz)
{   *ix = p->x[0];
    *iy = p->x[1];
    *iz = p->x[2];
}

inline unsigned getHashedIndex(unsigned ix, unsigned iy, unsigned iz)
{   return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK);
}

unsigned particle_bucket(particle_t* const p)
{
    unsigned ix, iy, iz;
    getIndices(p, &ix, &iy, &iz);

    return getHashedIndex(ix, iy, iz);
}

unsigned particle_neighborhood(unsigned* buckets, particle_t* const p)
{
    /* BEGIN TASK */

    // Get coordinates of the particl bine
    unsigned ix, iy, iz;
    getIndices(p, &ix, &iy, &iz);

    int counter = 0;
    int i,j,k;
    for(i = -1; i <= 1; i++)
    	{   
		int x = ix + i;
    		for(j = -1; j <= 1; j++)
    			{   
				int y = iy + j;
    				for(k = -1; k <= 1; k++)
  					{   
						int z = iz + k;
        					buckets[counter] = getHashedIndex(x,y,z);
        					counter++;
    					}
    			}	
    	}
    return counter;
    /* END TASK */
}

void hash_particles(sim_state_t* s)
{
    //HACK: Cleans hash, I think it can be done in a cleverer way without the mallocs
    free(s->hash);
    s->hash  = (particle_t**) malloc(HASH_SIZE*sizeof(particle_t*));
    int i;
    for (i = 0; i < HASH_SIZE; ++i)
        s->hash[i] = NULL;

    /* BEGIN TASK */
    particle_t* particle = s->part;
    while(particle != NULL) {
    	// Figure out bin for particle
        unsigned ix, iy, iz;
        getIndices(particle, &ix, &iy, &iz);
        int b = getHashedIndex(ix, iy, iz);

        // Add particle to the start of the list for bin b
        particle->next = s->hash[b];
        s->hash[b] = particle;

    	particle = particle->ll_next;
    }


    //TESTING
    /*
    for (i = 0; i < 10; ++i)
    {
        particle_t* part = &(s->part[i]);

        // Figure out bin for particle i
        unsigned ix, iy, iz;
        getIndices(part, &ix, &iy, &iz);
        int b = getHashedIndex(ix, iy, iz);
        printf("Hashindex for p%d (x:%u y:%u z:%u) is %d\n",i,ix, iy, iz, b);
        // Add particle to the start of the list for bin b
        part->next = s->hash[b];
        s->hash[b] = part;
    }
    */

}
