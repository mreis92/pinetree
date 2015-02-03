#include <stdio.h>
#include <string.h>

#include "align.h"

#define SEED(n,i,g,s) ((n-i-g >= s->seed_start && n-i-g <= s->seed_stop) ? 1 : 0)
#define CENTRAL_REGION(n,i,g,s) ((n-i-g >= s->cr_start && n-i-g <= s->cr_stop) ? 1 : 0)
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/* Read the file containing the results of FASTA execution */
fasta_info *process_alignment(char* filename){
	fasta_info *info;
	fasta_align **aligns;
	int miRNA_id, target_id;
	char miRNA_seq[BUFSIZE], target_seq[BUFSIZE];
	int start;
	int num_entries = 0;
	
	FILE *file = safe_fopen(filename, "r");
	int count = 0; 
	
	fscanf(file,"%d\n", &num_entries);
	
	info = (fasta_info*) safe_malloc(sizeof(fasta_info));
	aligns = (fasta_align**) safe_malloc(sizeof(fasta_align*)*num_entries);

	/* skip header */
	safe_free(get_string(file, BUFSIZE)); 
	
	while (fscanf(file,"%d,%d,%d,%[^,],%s\n", &miRNA_id,&target_id,&start,miRNA_seq, target_seq) != EOF){
		aligns[count] = (fasta_align*) safe_malloc(sizeof(fasta_align));
		
		aligns[count]->miRNA_id = miRNA_id;
		aligns[count]->target_id = target_id;
		aligns[count]->start = start;
		aligns[count]->target_seq = strdup(target_seq);
		aligns[count]->miRNA_seq = strdup(miRNA_seq);
		
		count++;
	}
	
	safe_fclose(file);
	
	info->aligns = aligns;
	info->count = count;
	
	return info;
}

/* Destroys the structure that contains the information from FASTA execution */
void clean_alignments(fasta_info* f){
	int i;
		
	for(i = 0; i < f->count; i++){
		safe_free(f->aligns[i]->target_seq);
		safe_free(f->aligns[i]->miRNA_seq);
		safe_free(f->aligns[i]);
	}
	
	safe_free(f->aligns);
	safe_free(f);
}

/* Attributes score to the alignment given by "align1" and "align2" */
float align_score(char* align1, char* align2, score_t *smodel){
	int i;
	uint n = strlen(align1);

	float match = 0;
	float mismatch = 0;
	float gap = 0;
	float gu = 0;
	
	float seed_mismatch = 0;
	float seed_gap = 0;
	float seed_gu = 0;
	
	/* shift seed one position for each gap found */
	int g = 0; 
	
	for(i = n-1; i >= 0; i--){
		float score;
		boolean seed;
		
		char a = align1[i];
		char b = align2[i];
		
		/* If gap in miRNA sequence, shift seed one nucleotide */
		if(b == '-')
			g++;
			
		seed = SEED(n,i,g,smodel);
		
		if((a == '-') || (b == '-')){
			if(seed)
				seed_gap++;
			else
				gap++;
			
			continue;
		}
		
		score = smodel->score_matrix[(int)a][(int)b];
		if(score == smodel->match)
			match++;
		else {
			if(score == smodel->wobble){
				if(seed)
					seed_gu++ ;
				else
					gu++;
			}
			else {
				if(score == smodel->mismatch){
					if(seed)
						seed_mismatch++;
					else
						mismatch++;
				}
			}
		}
	}
	
	return gap + mismatch + (seed_mismatch * 2) + (seed_gap * 2) + seed_gu + (gu * 0.5);
}

/* Outputs the kind of regulation mechanism for the miRNA:target pair */
char *mechanism(char* align1, char* align2, score_t *smodel){
	int i;
	uint n = strlen(align1);
	int center = smodel->cr_stop - smodel->cr_start + 1;

	/*FIXME: shift center one position for each gap found in transcript */
	int g = 0; 
	
	for(i = n-1; i >= 0; i--){
		boolean cr = CENTRAL_REGION(n,i,g,smodel);
		char a = align1[i];
		char b = align2[i];
		
		if(cr){
			float score = smodel->score_matrix[(int)a][(int)b];
			if((score != smodel->match) || (b == '-'))
				return "Translational inhibition";

			center--;
		}
		
		if(!center)
			return "Cleavage";
	}
	error("Unexpected error in mechanism detection");
}

/* Prints the alignment between align1 and align2 */
char* alignment_string(char* align1, char* align2, score_t *smodel){
	uint i;
	uint n = strlen(align1);
	char* alignment = (char*)safe_malloc(sizeof(char) * (n+1));

	for(i = 0; i < n; i++){
		float score;
		char a = align1[i];
		char b = align2[i];
		
		if((a == '-') || (b == '-')){
			alignment[i] = '-';
			continue;
		}
		
		score = smodel->score_matrix[(int)a][(int)b];
		if(score == smodel->match)
			alignment[i] = '|';
		else {
			if(score == smodel->wobble)
				alignment[i] = 'o';
			else {
				if(score == smodel->mismatch)
					alignment[i] = '.';
			}
		}
	}
	alignment[i] = '\0';
	return alignment;
}

/* Creates the structure that represents a scoring schema */
score_t *create_score_schema(float match, float mismatch, float gap, float wobble,
		     float seed_penalty, uint seed_start, uint seed_stop, uint cr_start, uint cr_stop)
{
	score_t *model = NULL;
	int i = 0, j = 0;
	float **score_matrix = NULL;

	model = (score_t *) safe_malloc(sizeof(score_t));
	model->gap = gap;
	model->seed_penalty = seed_penalty;
	model->match = match;
	model->mismatch = mismatch;
	model->wobble = wobble;
	model->seed_start = seed_start;
	model->seed_stop = seed_stop;
	model->cr_start = cr_start;
	model->cr_stop = cr_stop;

	score_matrix = (float **)safe_malloc(sizeof(float *) * 256);

	for (i = 0; i < 256; i++) {
		score_matrix[i] = (float *)safe_malloc(sizeof(float) * 256);
		for (j = 0; j < 256; j++) {
			score_matrix[i][j] = mismatch;
		}
	}

	score_matrix['a']['a'] = match;
	score_matrix['t']['t'] = match;
	score_matrix['A']['A'] = match;
	score_matrix['T']['T'] = match;
	score_matrix['u']['u'] = match;
	score_matrix['U']['U'] = match;
	
	score_matrix['c']['c'] = match;
	score_matrix['g']['g'] = match;
	score_matrix['G']['G'] = match;
	score_matrix['C']['C'] = match;

	score_matrix['G']['A'] = wobble;
	score_matrix['g']['a'] = wobble;
	score_matrix['T']['C'] = wobble;
	score_matrix['t']['c'] = wobble;
	score_matrix['U']['C'] = wobble;
	score_matrix['u']['c'] = wobble;

	model->score_matrix = score_matrix;

	return model;
}

/* Destroys the structure that represents a scoring schema */
void destroy_score_schema(score_t * model)
{
	int i = 0;

	if (model == NULL)
		return;

	for (i = 0; i < 256; i++)
		safe_free(model->score_matrix[i]);

	safe_free(model->score_matrix);
	safe_free(model);
}
