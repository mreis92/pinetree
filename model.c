#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "model.h"

stat_model_t *create_model(int markov_order){
	uint i;
	stat_model_t *model = (stat_model_t *)safe_malloc(sizeof(stat_model_t));
	model->markov_order = markov_order;
	model->entries = pow(NUCLEOTIDES,markov_order)+1;
	model->p = (float*)safe_calloc(model->entries, sizeof(float));
	model->A = (float**)safe_malloc(sizeof(float*) * model->entries);
	
	for(i = 0; i < model->entries; i++){
		model->A[i] = (float*)safe_calloc(NUCLEOTIDES+1, sizeof(float));
	}
	
	return model;
}

void destroy_model(stat_model_t *model){
	uint i;
	
	safe_free(model->p);
	
	for(i = 0; i < model->entries; i++){
		safe_free(model->A[i]);
	}
	
	safe_free(model->A);
	safe_free(model);
}

/* Creates a set of models from a dataset. Each sequence
    of the dataset creates a new statistical model */
stat_models_t *create_gene_models(dataset_t * dataset, int markov_order){
	uint i;
	stat_models_t *models = (stat_models_t *)safe_malloc(sizeof(stat_models_t));

	models->modeln = dataset->seqn;
	models->models = (stat_model_t **) safe_malloc(sizeof(stat_model_t *) * models->modeln);
	models->ids = dataset->ids;

	for (i = 0; i < models->modeln; i++) {
		models->models[i] = create_model(markov_order);
		update_model(models->models[i], dataset->sequences[i]);
		normalize_model(models->models[i]);
	}

	return models;
}

/* Destroys the set of models given by "models" */
void destroy_gene_models(stat_models_t * models){
	uint i;

	for (i = 0; i < models->modeln; i++)
		destroy_model(models->models[i]);

	safe_free(models->models);
	safe_free(models);
}

/* Creates a model of the background, from the set of sequences given
    by "background_ds". Unlike a gene model, each of the sequences in
    the dataset will contribute to the same model               */
stat_model_t *create_background_model(dataset_t * background_ds, int markov_order){

	stat_model_t *model = create_model(markov_order);
	uint i;
	for (i = 0; i < background_ds->seqn; i++) {
		update_model(model, background_ds->sequences[i]);
	}
	normalize_model(model);

	return model;
}

/* Given a model and a sequence, it updates the model with the
    new sequence. This ways, there is no transition between the
    last character of a sequence and the first character of the
    next sequence                                               */
/*TODO: Refactor repeated code. */
void update_model(stat_model_t * model, char *seq){
	int i;
	int stop = model->markov_order;
	int line_length = strlen(seq);
	int count_index = model->entries-1;
	long long int mask = 0;
	long long int previous = 0;
	
	for(i = 0; i < stop; i++){
		int temp = bin_code(seq[i]);
		if(temp != 0){
			mask = (mask<<2)|0b11;
			previous = ((previous<<2)&mask) |(temp&0b011);
		}
		else {
			mask = 0;
			stop = model->markov_order+i+1;
		}
	}
 
	model->p[previous]++;
	model->p[count_index]++;

	for (; i < line_length; i++) {
		int cur_char = bin_code(seq[i]);
		llint pos;
		
		if (cur_char != 0) {
			cur_char = cur_char&0b011;
			pos = ((previous<<2)&mask)|cur_char;
			model->p[pos]++;
			model->p[count_index]++;
		
			model->A[previous][cur_char]++;
			model->A[previous][_N_]++;
			
			previous = pos;
			
		}
		else{
			stop = model->markov_order+(++i);
			
			for(; (i < stop) && (i < line_length); i++){
				cur_char = bin_code(seq[i]);
				if(cur_char != 0){
					previous = ((previous<<2)&mask)|(cur_char&0b011);
				}
				else 
					stop = model->markov_order+i+1;
			}
			
			if(i < line_length){
				model->p[previous]++;
				model->p[count_index]++;
				i--;
			}
		}
	}
}

/* Aplies pseudocounts to the model and changes the
    representation from frequencies to probabilities.   */
void normalize_model(stat_model_t * model){
	int i, j;
	int count_index = model->entries-1;

	pseudocounts(model);

	for (i = 0; i < count_index; i++) {
		model->p[i] /= model->p[count_index];
	}

	
	for (i = 0; i < count_index; i++) {
		for(j = 0; j < NUCLEOTIDES; j++){
			model->A[i][j] /= model->A[i][_N_];
		}
	}
}

/* Applies pseudocounts to the Markov chains of a statistical model */
void pseudocounts(stat_model_t * model)
{
	int i, j;
	int count_index = model->entries-1;
	
	for (i = 0; i < count_index; i++) {
		model->p[i] += PSEUDOCOUNTS;
		model->p[count_index] += PSEUDOCOUNTS;
	}
	
	for (i = 0; i < count_index; i++) {
		for(j = 0; j < NUCLEOTIDES; j++){
			model->A[i][_N_] += PSEUDOCOUNTS;
			model->A[i][j] += PSEUDOCOUNTS;
		}
	}
}

/* Prints relevant information about a statistical model */
void print_model(stat_model_t * model){
	int i, j;
	int count_index = model->entries-1; 
	char *previous = (char*)safe_malloc((sizeof(char)*model->markov_order)+1);

	for (i = 0; i < count_index; i++) {
		llint temp = i;
				
		for(j = (model->markov_order-1); j >= 0; j--){
			previous[j] = char_code(temp&0b11);
			temp >>=2;
		}
		previous[model->markov_order] = '\0';
			
		printf("p(%s) = %f\n", previous, model->p[i]);
	}

	for (i = 0; i < count_index; i++) {
		llint temp = i;
		
		for(j = (model->markov_order-1); j >= 0; j--){
			previous[j] = char_code(temp&0b11);
			temp >>=2;
		}
		previous[model->markov_order] = '\0';
		
		for (j = 0; j < NUCLEOTIDES; j++) {
			printf("p(%s)(%c) = %f\n", previous, char_code(j), model->A[i][j]);
		}
	}
	printf("\n");
	
	safe_free(previous);
}

/* Initializes an entity that contains information about a mirna in this 
	context. Each mirna, will have a score for each gene and background, with a certain number
	of allowed errors */
mirna_info_t *create_mirna_info(uint gene_seqn)
{
	int e;

	mirna_info_t *mirna = (mirna_info_t *)safe_malloc(sizeof(mirna_info_t));
	mirna->gene_prob =
	    (float **)safe_malloc(sizeof(float *) * (NUM_ERRORS+1));
	mirna->background_prob =
	    (float *)safe_malloc(sizeof(float) * (NUM_ERRORS+1));

	for (e = 0; e <= NUM_ERRORS; e++)
		mirna->gene_prob[e] =
		    (float *)safe_malloc(sizeof(float) * gene_seqn);

	return mirna;
}

mirnas_info_t *create_mirnas_info(uint seqn, uint gene_seqn)
{
	uint i;
	mirnas_info_t *mirnas = (mirnas_info_t *)safe_malloc(sizeof(mirnas_info_t));
	mirnas->mirna = (mirna_info_t **)safe_malloc(sizeof(mirna_info_t *) * seqn);
	mirnas->seqn = seqn;

	for (i = 0; i < seqn; i++) {
		mirnas->mirna[i] = create_mirna_info(gene_seqn);
	}

	return mirnas;
}

void destroy_mirna_info(mirna_info_t * mirna)
{
	int e;
	for (e = 0; e <= NUM_ERRORS; e++)
		safe_free(mirna->gene_prob[e]);

	safe_free(mirna->gene_prob);
	safe_free(mirna->background_prob);
	safe_free(mirna->sequence); /* it's the reverse complementar */
	safe_free(mirna);
}

void destroy_mirnas_info(mirnas_info_t * mirnas)
{
	uint i;

	for (i = 0; i < mirnas->seqn; i++) {
		destroy_mirna_info(mirnas->mirna[i]);
	}

	safe_free(mirnas->mirna);
	safe_free(mirnas);
}
