#ifndef EVO_H
#define EVO_H

#include "types.h"
#include "util.h"
#include "dataset.h"
#include "model.h"
#include "fasta.h"
#include "constants.h"

typedef struct evo_info_t {
	stat_model_t *background_model;
	stat_models_t *gene_models;
	mirnas_info_t *mirnas;
} evo_info_t;

ldouble monad_errors(stat_model_t * model, char *seq, uint num_errors, int offset, int limit, llint previous, llint mask);
ldouble prob_monad_errors(stat_model_t * model, char *seq, uint num_errors);
void prob_score(stat_model_t * genome_model,
		stat_models_t * gene_models, mirnas_info_t * mirnas, int np);
ldouble calculate_escore(evo_info_t* evo_info, int mirna_id, int gene_id);
evo_info_t* generate_models(dataset_t *background_seqs, dataset_t *gene_seqs, dataset_t *mirna_seqs, int markov_order, int np);
void cleanup(evo_info_t* evo_info);

#endif
