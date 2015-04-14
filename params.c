#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "params.h"

static int* set_n_et(sim_param_t* params);
static float* randomize_dconst(sim_param_t* params);
static float* get_rate_constants(sim_param_t* params);
static void get_diffusion_constants(sim_param_t* params);

//hardcoded some default parameters for testing
void default_params(sim_param_t* params) {
	params->fname = "out.txt";
	params->N = 5;
	params->M_uni = 1;
	params->M_bi = 0;
	params->n_t = 2;
	params->n_et = set_n_et(params);
	params->dim_size = 100.0;
	params->cap_r = 50.0;
	//params->dconst = randomize_dconst(params);
	params->rateconstants = get_rate_constants(params);
	//params->diff_dist_h = 5;
	get_diffusion_constants(params);
}

static void get_diffusion_constants(sim_param_t* params) {
	float* constants = (float*) calloc(params->n_t, sizeof(float));
	float* diffdistconst = (float*) calloc(params->n_t, sizeof(float));

	FILE* file = fopen("diffusion.txt","r");
	char line[256];
	int lineno = 0;
	int flag = 0;//flag to tell whether we are reading diffusion constants or diffusion distant constants
	int dconstcount = 0;
	int diffdistcount = 0;

	while(fgets(line, sizeof(line), file)) {
		if(!flag) {
			if(dconstcount == params->n_t ) {
				flag = 1;
			} else {
				double value = atof(line);
				constants[dconstcount] = value;
				dconstcount++;
			}
		}
		if(flag) {
			double value = atof(line);
			diffdistconst[diffdistcount] = value;
			diffdistcount++;

		}
	}

	params->dconst = constants;
	params->diff_dist_h = diffdistconst;
}

static float* get_rate_constants(sim_param_t* params) {
	float* constants = (float*) calloc(params->M_uni + params->M_bi, sizeof(float));

	//NOTE:in the file, all unimolecular reactions' constants must come first, then bi
	FILE* file = fopen("rateconstants.txt", "r");
	char line[256];
	int lineno = 0;

	while(fgets(line, sizeof(line), file)) {
		double value = atof(line);
		constants[lineno] = (float) value;

		lineno++;
	}
	/*
	int i;
	printf("rateconstants:\n");
	for(i=0;i<params->M_uni + params->M_bi;++i) {
		printf("%f\n",constants[i]);
	}
	*/
	return constants;
}

static float* randomize_dconst(sim_param_t* params) {
	float* dconst = (float*) calloc(params->n_t, sizeof(float));

	int i;
	for(i=0;i<params->n_t;++i) {
		dconst[i] = (double)rand() / (double)RAND_MAX;//random value between 0-1
	}
	return dconst;
}

static int* set_n_et(sim_param_t* params) {
	srand(time(0));

	int* a = (int*) calloc(params->n_t, sizeof(int));
	FILE* file = fopen("particlenos.txt","r");
	char line[256];
	int lineno=0,i;

	while(fgets(line, sizeof(line),file)) {
		a[lineno] = atoi(line);
		lineno++;
	}

	for(i=0; i<params->n_t;++i) {
		printf("numof(type%d):%d\n",i, a[i] );
	}

	return a;
	// MANUAL SETTING OF n_et
	//int* a = (int*) calloc(params->n_t, sizeof(int));
}

void get_reactions(sim_param_t* params) {
	//allocate space for reaction list
	params->reactions = (reaction*) calloc((params->M_uni + params->M_bi), sizeof(reaction));

	int max_input_reactions = 2;
	//parse the file and fill the reaction struct
	FILE* file = fopen("reactions.txt", "r");
	char line[256];
	int lineno = 0;

	while(fgets(line, sizeof(line), file)) {

		int i;

		int outp = 0, coeff = 1, pair = 0;

		int* inpt = (int*) malloc(max_input_reactions*sizeof(int));//temp array for input tyeps
		int* inpc = (int*) malloc(max_input_reactions*sizeof(int));//temp array for input coeff
		int inputindex = 0;
		int opt;//output type
		int opc;//output coeffecient
		int n=0, t=0; //coefficient and type pair

		for(i = 0; i < strlen(line); ++i) {
			int t1;
			if(line[i] == '\n') {
				break;
			}
			if (line[i] == '+') {
				continue;
			}
			if (line[i] == '=') {
				outp = 1;
				continue;
			}

			if(outp) {
				if(line[0] == '0') { //special case for unimolecular zero output condition
					opc = 0;
					break;
				}
				if((int)line[i] < 97) {
					n*=10;
					n += line[i] - '0';
				} else {
					coeff = 0;
				}
				if(!coeff) {
					t1 = (int) line[i];
					t = t1 - 97;
					opt = t;
					opc = n;
					coeff = 1;
				}
			} else {
				if((int)line[i] < 97) {
					n*=10;
					n += line[i] - '0';
				} else {
					coeff = 0;
				}
				if(!coeff) {
					t1 = (int) line[i];
					t = t1 - 97; // lower case a starts from int 97
					inpc[inputindex] = n;
					inpt[inputindex] = t;
					coeff = 1;
					inputindex++;
					n=0;
				}
			}
		}

		//allocate the reaction we have just parsed
		reaction* react = (reaction*) calloc(1, sizeof(reaction));
		react->input = (reactant*) calloc(inputindex, sizeof(reactant));

		//put input in react
		for(i = 0; i < inputindex; ++i) {
			react->input[i].n = inpc[i];
			react->input[i].ptype = inpt[i];
		}
		//put output in react
		react->output.n = opc;
		react->output.ptype = opt;

		//set whether reaction is uni or bimolecular
		if(inputindex == 1 ) {
			react->type = 1;
		} else {
			react->type = 2;
		}

		//finally add reaction(react) to reaction list
		params->reactions[lineno] = *react;


		free(inpt);
		free(inpc);
		lineno++;
	}

	fclose(file);
}

