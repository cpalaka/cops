#ifndef PARAMS_H
#define PARAMS_H


//holds info on reaction of type xA + yB -> zC where x,y,z are coefficients
//and A,B,C are particle types (xA, yB, and zC are held as reactants)
typedef struct reactant {
	int n;	  //coefficient
	int ptype;//particle type
} reactant;

typedef struct reaction {
	reactant* input;
	reactant output;
	int type; //1 - unimolecular, 2 - bimolecular
} reaction;

typedef struct sim_param_t {
	char* fname;    //output file name
	int N;          //number of particles
	int M_uni;	    //number of unimolecular reactions
	int M_bi;	    //number of bimolecular reactions
	int n_t;        //number of particle types
	int* n_et;      //number of particles of each type ( 0 - (n_t-1))
	reaction* reactions;//list of reactions
	float dim_size;   //dimension of sim area
	float cap_r;      //size of reaction capture radius
	float* rateconstants;//rate constants (uni comes first, then bi)
	float* dconst;     //diffusion constant
	float* diff_dist_h;  //constant distant to move particle in diffusion event
} sim_param_t;


void default_params(sim_param_t* params);
void get_reactions(sim_param_t* params);


#endif 