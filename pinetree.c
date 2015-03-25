#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <omp.h>
#include <pthread.h>
#include <time.h>
#include <sys/wait.h>
#include <unistd.h>

#include "accessibility.h"
#include "align.h"
#include "annotation.h"
#include "dataset.h"
#include "constants.h"
#include "fasta.h"
#include "pinetree_cmdl.h"
#include "pinetree_utils.h"
#include "strmap.h"
#include "util.h"

void print_pair_info(fasta_align *cur_align, pinetree_args *args, dataset_t *tds, dataset_t *mds, FILE* output_file){
	char *seq1 = cur_align->target_seq;
	char *seq2 = cur_align->miRNA_seq;
	int tid = cur_align->target_id;
	int mid = cur_align->miRNA_id;

	char* align = alignment_string(seq1, seq2, args->sschema);
	char* reg_mechanism = mechanism(seq1, seq2, args->sschema);
	
	if(args->human_output){
		fprintf(output_file, "target id: %s\n", tds->ids[tid]);

		if(tds->annotations)
			fprintf(output_file, "target info: %s\n", tds->annotations[tid] ? tds->annotations[tid] : "N/A");
			
		fprintf(output_file, "miRNA id: %s\n", mds->ids[mid]);
		fprintf(output_file, "Reg. mechanism: %s\n", reg_mechanism);

		fprintf(output_file, "Complementarity score: %.*f\n", 
			(args->normalization ? 5 : 1), cur_align->cscore);

		if(args->accessibility)
			fprintf(output_file, "Accessibility score: %.*f\n", 
				(args->normalization ? 5 : 1), cur_align->ascore);

		fprintf(output_file, "miRNA:\t%s\n", seq2);
		fprintf(output_file, "\t%s\n", align);
		fprintf(output_file, "target:\t%s\n", seq1);
		fprintf(output_file, "#\n");
	}	
	else {
		fprintf(output_file, "%s,%s,%.*f", tds->ids[tid], mds->ids[mid], 
			(args->normalization ? 5 : 1), cur_align->cscore);

		if(args->accessibility)
			fprintf(output_file, ",%.*f", 
				(args->normalization ? 5 : 1), cur_align->ascore);

		fprintf(output_file,",%s,%s", align, reg_mechanism);
	
		if(tds->annotations)
			fprintf(output_file, ",%s", tds->annotations[tid] ? tds->annotations[tid] : "N/A");
		fprintf(output_file, "\n");
	}
	
	safe_free(align);
}

void target_prediction(fasta_info *info, pinetree_args* args, dataset_t *tds, dataset_t *mds){
	
	uint i;
	fasta_align **aligns = info->aligns;
	
	#pragma omp parallel num_threads(args->num_processors)
	{
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");
		
		#pragma omp for 
		for(i = 0; i < info->count; i++){
			char *align1 = aligns[i]->target_seq;
			char *align2 = aligns[i]->miRNA_seq;
			aligns[i]->cscore = align_score(align1, align2, args->sschema);
			
			if(aligns[i]->cscore <= args->c_threshold){
				if(args->accessibility && (aligns[i]->ascore > args->a_threshold))
					continue;

				print_pair_info(aligns[i], args, tds, mds, output_file);
			}
			
		}
		
		safe_fclose(output_file);
	}
}

void target_prediction_norm(fasta_info *info, pinetree_args* args, dataset_t *tds, dataset_t *mds){
	uint i;
	fasta_align **aligns = info->aligns;
	float c_avg = 0.0, c_std_dev = 0.0;
	float a_avg = 0.0, a_std_dev = 0.0;
	
	#pragma omp parallel for reduction(+:c_avg, a_avg) num_threads(args->num_processors) 
	for(i = 0; i < info->count; i++){
		char *align1 = aligns[i]->target_seq;
		char *align2 = aligns[i]->miRNA_seq;
		float cscore = align_score(align1, align2, args->sschema);
		aligns[i]->cscore = cscore;
		c_avg += cscore;
		a_avg += aligns[i]->ascore;
	}

	c_avg /= info->count;
	a_avg /= info->count;

	#pragma omp parallel for reduction(+:c_std_dev, a_std_dev) num_threads(args->num_processors) 
	for(i = 0; i < info->count; i++){
		c_std_dev += (aligns[i]->cscore-c_avg) * (aligns[i]->cscore-c_avg);
		a_std_dev += (aligns[i]->ascore-a_avg) * (aligns[i]->ascore-a_avg);
	}

	c_std_dev = sqrt(c_std_dev/info->count);
	a_std_dev = sqrt(a_std_dev/info->count);

	#pragma omp parallel num_threads(args->num_processors)
	{
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");
		
		#pragma omp for 
		for(i = 0; i < info->count; i++){
			aligns[i]->cscore = (aligns[i]->cscore - c_avg)/c_std_dev;
			
			if(aligns[i]->cscore <= args->c_threshold){
				if(args->accessibility){
					aligns[i]->ascore = (aligns[i]->ascore - a_avg)/a_std_dev;
					
					if(aligns[i]->ascore > args->a_threshold)
						continue;
				}

				print_pair_info(aligns[i], args, tds, mds, output_file);
			}
		}
		
		safe_fclose(output_file);
	}
}

void calc_acc_values(pinetree_args *args, fasta_align *cur_align, uint current, dataset_t *tds, dataset_t *mds, FILE* temp_file){
	char *seq1 = cur_align->target_seq;
	char *seq2 = cur_align->miRNA_seq;
	int tid = cur_align->target_id;
	int mid = cur_align->miRNA_id;
	int start = cur_align->start;
 	int length = strlen(mds->sequences[mid]) - count_occurrences(seq1, '-');
	
	if(!args->normalization){
		float cscore = align_score(seq1, seq2, args->sschema);

		if(cscore > args->c_threshold)
			return;
	}

	double ascore = calculate_accessibility("RNAup", tds->sequences[tid], start, length);

	/* FIXME */
	if(ascore == INFINITY)
		ascore = 25;

	fprintf(temp_file, "%d,%.1f\n", current, ascore);
}

void load_acc_values(pinetree_args *args, fasta_align **aligns){

	#pragma omp parallel num_threads(args->num_processors)
	{
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "r");

		uint current = 0;
		float ascore = 0.0;

		while (fscanf(output_file,"%u,%f\n", &current, &ascore) != EOF)
			aligns[current]->ascore = ascore;

		safe_fclose(output_file);
		safe_remove(args->temp_file[omp_get_thread_num()]);
	}

}

int main(int argc, char **argv){
	uint j;
	
	pinetree_args* args = read_cml_arguments(argc, argv);
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	
	StrMap *target_map = map_ids(tds);
	StrMap *miRNA_map = map_ids(mds);
	
	fasta_info *info = process_alignment(args->num_processors, args->evalue, 
									args->mirna_file, args->transcript_file, miRNA_map, target_map);
	
	if(args->annotation_file)
		annotate_targets(tds, target_map, args->annotation_file);
	
	sm_delete(target_map);
	sm_delete(miRNA_map);

	/*We have to use this ugly fork method instead of using OpenMP beacuse RNAup doesn't 
	 * like when we use parallellize code in the outer execution... probably some
	 static variables exist, and to have serial analysis of accessibility is not feasible.

	 Right now accessibility calculation is a two-step process -- write values to disk, and
	 then read them to memory. This makes it easier to perform the parallel normalization
	 if needed. */
	if(args->accessibility){
		uint id;
		int status;
		pid_t wpid;
		int batch = info->count/args->num_processors;
		int remainder = info->count%args->num_processors;

		for (id = 0; id < args->num_processors; ++id) {
			if ((wpid = fork()) < 0) 
				printf("Error in fork\n");
			else if (wpid == 0) {
				uint i;
				uint lower = id*batch;
				uint upper = (id==args->num_processors-1) ? (id+1)*batch + remainder : (id+1)*batch;
				FILE *temp_file = safe_fopen(args->temp_file[id], "w");
				
				for (i = lower; i < upper; i++) 
					calc_acc_values(args, info->aligns[i], i, tds, mds, temp_file);
				
				safe_fclose(temp_file);
				exit(0);
			}
		}
		while ((wpid = wait(&status)) > 0); /* Wait for children to exit. */
		load_acc_values(args, info->aligns);
	}

	if(args->normalization)
		target_prediction_norm(info, args, tds, mds);
	else /* default mode */
		target_prediction(info, args, tds, mds);
	
	file_joiner(args);

	clean_alignments(info);
	destroy_dataset(tds);
	destroy_dataset(mds);
	destroy_score_schema(args->sschema);
	
	safe_free(args);

	return 0;
}



