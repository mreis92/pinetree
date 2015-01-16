#ifndef ALIGN_H
#define ALIGN_H

#include "types.h"
#include "util.h"

typedef struct score_t {
	float match;
	float wobble;
	float mismatch;
	float gap;		/* gap penalty */
	uint cr_start;
	uint cr_stop;
	uint seed_start;
	uint seed_stop;
	float seed_penalty;	/* seed penalty */
	float **score_matrix;
} score_t;

float align_score(char* align1, char* align2, score_t *smodel);
char* mechanism(char* align1, char* align2, score_t *smodel);
char* alignment_string(char* align1, char* align2, score_t *smodel);

score_t *create_score_schema(float, float, float, float, float, uint, uint, uint, uint);
void destroy_score_schema(score_t *);

#endif
