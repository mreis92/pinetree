#ifndef ALIGN_H
#define ALIGN_H

#include <stdio.h>

#include "types.h"
#include "util.h"

typedef struct fasta_align{
	char* miRNA_seq;
	char* target_seq;
	int target_id;
	int miRNA_id;
	int start;
	float cscore;
	float ascore;
} fasta_align;

typedef struct fasta_info{
	fasta_align **aligns;
	int count;
} fasta_info;

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

fasta_info *process_alignment(char* filename);
void clean_alignments(fasta_info* f);
float align_score(char* align1, char* align2, score_t *smodel);
char* mechanism(char* align1, char* align2, score_t *smodel);
char* alignment_string(char* align1, char* align2, score_t *smodel);

score_t *create_score_schema(float, float, float, float, float, uint, uint, uint, uint);
void destroy_score_schema(score_t *);

#endif
