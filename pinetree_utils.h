#ifndef PINETREE_UTILS_H
#define PINETREE_UTILS_H

#include <stdio.h>

#include "align.h"
#include "constants.h"
#include "types.h"
#include "util.h"

typedef struct pinetree_args {
	score_t *sschema;
	uint num_processors;
	char* transcript_file;
	char* mirna_file; 
	char* align_file; 
	char* annotation_file;
	char* output_file;
	char** temp_file;
	char* start_time;
	char param_info[LONGBUF];
	float c_threshold;
	float a_threshold;
	float e_threshold;
	uint num_errors;
	uint markov_order;
	boolean normalization;
	boolean human_output;
	boolean accessibility;
} pinetree_args;


void file_joiner(pinetree_args* args, char* header);

#endif
