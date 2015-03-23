#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "evo.h"

/* Recursive function for calculation of the probability appearing 
	in the Markov chain of a statistical model */
float monad_errors(stat_model_t * model, char *seq, uint num_errors, int offset, int limit, llint previous, llint mask){
	int i;
	float probability = 0;
	int temp = (bin_code(seq[offset]))&0b011;
	
	if (offset < limit) {
		for (i = 0; i < NUCLEOTIDES; i++){
			llint current = ((previous<<2)&mask)|i;
			
			if (i != temp && num_errors > 0){
				probability += model->A[previous][i] * monad_errors(model, seq, num_errors - 1, offset + 1, limit, current, mask);
			}
			else if (i == temp) {
				probability += model->A[previous][i] * monad_errors(model, seq, num_errors, offset + 1, limit, current, mask);
			}
		}
	} else {
		for (i = 0; i < NUCLEOTIDES; i++)
			if ((i != temp && num_errors > 0) || i == temp) {
				probability += model->A[previous][i];
			}
	}

	return probability;
}

uint count_errors(llint a, llint b, int order, uint num_errors){
	int i;
	uint errors = 0;
	
	for(i = 0; i < order; i++){
		if((a&0b11)!=(b&0b11)){
			errors++;
		
			if(errors > num_errors)
				return errors;
		}
			
		a>>=2;
		b>>=2;
	}
	
	return errors;
	
}
/* Calculates the probability of a sequence appearing in the Markov chain given by "model", up to "num_errors" */
float prob_monad_errors(stat_model_t * model, char *seq, uint num_errors){
	llint i, current = 0, mask = 0;
	float probability = 0;
	int limit = strlen(seq) - 1;
	int entries = model->entries-1;
	int stop = model->markov_order;
	
	
	for(i = 0; i < stop; i++){
		int temp = bin_code(seq[i]);
		if(temp != 0){
			mask = (mask<<2)|0b11;
			current = ((current<<2)&mask) |(temp&0b011);
		}
		else {
			mask = 0;
			stop = model->markov_order+i+1;
		}
	}

	for (i = 0; i < entries; i++) {
		uint errors = count_errors(i, current, model->markov_order, num_errors); 
		
		if((errors > 0) && (errors <= num_errors))
			probability += model->p[i] * monad_errors(model, seq, num_errors-errors, model->markov_order, limit, i, mask);
		else if(errors == 0){
			probability += model->p[i] * monad_errors(model, seq, num_errors, model->markov_order, limit, i, mask);
		}
	}

	return probability;
}

/* Calculates the probability of a miRNA appearing in the background and genes,
    according to their Markov models                                    */
void prob_score(stat_model_t * background_model, stat_models_t * gene_models, mirnas_info_t * mirnas, int np){
	uint i, j, e;
	omp_set_num_threads(np);

	for (e = 0; e <= NUM_ERRORS; e++) {
		#pragma omp parallel for private(j)
		for (i = 0; i < mirnas->seqn; i++) {
			mirnas->mirna[i]->background_prob[e] = prob_monad_errors(background_model, mirnas->mirna[i]->sequence, e);

			for (j = 0; j < gene_models->modeln; j++) {
				mirnas->mirna[i]->gene_prob[e][j] = prob_monad_errors(gene_models->models[j], mirnas->mirna[i]->sequence, e);
				
			}
		}
	}
}

/* Prints the scores for each miRNA against the background and genes */
float calculate_affinity(evo_info_t* evo_info, int mirna_id, int gene_id)
{
	float background_prob, gene_prob;

	background_prob = evo_info->mirnas->mirna[mirna_id]->background_prob[0];
	gene_prob = evo_info->mirnas->mirna[mirna_id]->gene_prob[0][gene_id];

	return logl(gene_prob) - logl(background_prob);
}

evo_info_t* generate_models(dataset_t *background_seqs, dataset_t *gene_seqs, dataset_t *mirna_seqs, int markov_order, int np)
{
	uint i;

	evo_info_t* evo_info = (evo_info_t*) malloc(sizeof(evo_info_t));

	/* Create models from datasets */
	evo_info->background_model = create_background_model(background_seqs, markov_order);
	evo_info->gene_models = create_gene_models(gene_seqs, markov_order);

	evo_info->mirnas =
	    create_mirnas_info(mirna_seqs->seqn, evo_info->gene_models->modeln);

	/* Reverse mirnas for logit calculation */
	for (i = 0; i < mirna_seqs->seqn; i++) {
		evo_info->mirnas->mirna[i]->sequence = reverse_complement(mirna_seqs->sequences[i]);
		evo_info->mirnas->mirna[i]->id = mirna_seqs->ids[i];
		evo_info->mirnas->mirna[i]->length = strlen(evo_info->mirnas->mirna[i]->sequence); /*TODO: inneficient */
	}

	/* Calculate prob and print scores */
	prob_score(evo_info->background_model, evo_info->gene_models, evo_info->mirnas, np);
	/*print_logit_scores(mirnas, background_model, gene_models);*/


	return evo_info;
}

void cleanup(evo_info_t* evo_info){
	safe_free(evo_info->background_model);
	destroy_gene_models(evo_info->gene_models);
	destroy_mirnas_info(evo_info->mirnas);
	safe_free(evo_info);
}

