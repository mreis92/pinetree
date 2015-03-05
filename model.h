#ifndef MODEL_H
#define MODEL_H

#include "constants.h"
#include "dataset.h"
#include "fasta.h"
#include "types.h"
#include "util.h"

/* A statistical model that contain the Markov chains of order 0 and 1 for a given sequence */
typedef struct stat_model_t {
	uint markov_order;
	uint entries;
	ldouble *p;
	ldouble **A;
} stat_model_t;

/* A set of statistical models, with an id associated with each model */
typedef struct stat_models_t {
	uint modeln;
	stat_model_t **models;
	char **ids;
} stat_models_t;

/* The relevant information about a mirna, namely its sequence, id, and its background
	probability + probability of appearing on a given transcript */
typedef struct mirna_info_t {
	char *sequence;
	char *id;
	uint length;
	ldouble *background_prob;
	ldouble **gene_prob;

} mirna_info_t;

/* Information about a set of mirnas */
typedef struct mirnas_info_t {
	mirna_info_t **mirna;
	uint seqn;
} mirnas_info_t;

stat_model_t *create_model(int markov_order);
void destroy_model(stat_model_t *model);
stat_models_t *create_gene_models(dataset_t * dataset, int markov_order);
void destroy_gene_models(stat_models_t * models);
stat_model_t *create_background_model(dataset_t * background_ds, int markov_order);
void update_model(stat_model_t * model, char *seq);
void normalize_model(stat_model_t * model);
void pseudocounts(stat_model_t * model);
void print_model(stat_model_t * model);

mirna_info_t *create_mirna_info(uint gene_seqn);
mirnas_info_t *create_mirnas_info(uint seqn, uint gene_seqn);
void destroy_mirna_info(mirna_info_t * mirna);
void destroy_mirnas_info(mirnas_info_t * mirnas);

#endif
