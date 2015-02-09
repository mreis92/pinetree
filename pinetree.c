#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <omp.h>
#include <time.h>
#include <sys/wait.h>
#include <unistd.h>

#include "accessibility.h"
#include "align.h"
#include "annotation.h"
#include "dataset.h"
#include "constants.h"
#include "fasta.h"
#include "pinetree_utils.h"
#include "util.h"

/*Fairly big scores to initialize variables for normalization. */
#define MAX_CSCORE 50.0
#define MIN_CSCORE 0.0

pinetree_args *initialize_args(){
	pinetree_args *args = (pinetree_args*)safe_malloc(sizeof(pinetree_args));
	
	args->sschema = create_score_schema(MATCH, MISMATCH, GAP, GU, SEED_PENALTY, SEED_START, SEED_STOP, CR_START, CR_STOP);
	args->num_processors = NUM_PROCESSORS;
	args->transcript_file = NULL;
	args->mirna_file = NULL;
	args->align_file = NULL;
	args->output_file = NULL;
	args->annotation_file = NULL;
	args->c_threshold = C_THRESHOLD;
	args->a_threshold = A_THRESHOLD;
	args->start_time = 0;
	args->accessibility = 0;
	args->normalization = 0;
	args->human_output = 0;
	
	return args;
}

void initialize_parameters(char* filename, pinetree_args* args){
	char buffer[BUFSIZE];
	char *parameter, *value;

	FILE *config_file = safe_fopen(filename, "r");
	score_t *s = args->sschema;

	while(fgets (buffer, BUFSIZE, config_file) != NULL){
		parameter = strtok(buffer, "=");
		value = strtok(NULL, "");

		if(STRMATCH("NUM_PROC", parameter))
			sscanf(value, "%u", &(args->num_processors));
		else if(STRMATCH("GAP", parameter))
			sscanf(value, "%f", &(s->gap));
		else if(STRMATCH("MATCH", parameter))
			sscanf(value, "%f", &(s->match));
		else if(STRMATCH("MISMATCH", parameter))
			sscanf(value, "%f", &(s->mismatch));
		else if(STRMATCH("GU_WOBBLE", parameter))
			sscanf(value, "%f", &(s->wobble));
		else if(STRMATCH("SEED_PENALTY", parameter))
			sscanf(value, "%f", &(s->seed_penalty));
		else if(STRMATCH("SEED_REGION_LOWER_BOUND", parameter))
			sscanf(value, "%u", &(s->seed_start));
		else if(STRMATCH("SEED_REGION_UPPER_BOUND", parameter))
			sscanf(value, "%u", &(s->seed_stop));
		else if(STRMATCH("CENTRAL_REGION_LOWER_BOUND", parameter))
			sscanf(value, "%u", &(s->cr_start));
		else if(STRMATCH("CENTRAL_REGION_UPPER_BOUND", parameter))
			sscanf(value, "%u", &(s->cr_stop));
		else if(STRMATCH("COMP_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->c_threshold));
		else if(STRMATCH("ACC_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->a_threshold));
		else if(STRMATCH("NCOMP_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->nc_threshold));
		else if(STRMATCH("NACC_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->na_threshold));
	}
	
	safe_fclose(config_file);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	char buffer[BUFSIZE];
	int i, c;
	int tflag = 0, mflag = 0, aflag = 0;
	time_t mytime = time(NULL);
	char *config_file = NULL;
	pinetree_args *args = initialize_args();
	score_t *s = args->sschema;

	opterr = 0;
	
	while ((c = getopt (argc, argv, "C:t:m:a:o:n:A:s:zNH")) != -1){
		switch (c){
		case 'C':
			config_file = optarg;
			initialize_parameters(config_file, args);
			break;
		case 't':
			tflag = 1;
			args->transcript_file = optarg;
			break;
		case 'm':
			mflag = 1;
			args->mirna_file = optarg;
			break;
		case 'a':
			aflag = 1;
			args->align_file = optarg;
			break;
		case 'o':
			args->output_file = optarg;
			break;
		case 'n':
			args->num_processors = atoi(optarg);
			break;
		case 'A':
			args->annotation_file = optarg;
			break;
		case 's':
			mytime = atoi(optarg);
			break;
		case 'z':
			args->accessibility = 1;
			break;
		case 'N':
			args->normalization = 1;
			break;
		case 'H':
			args->human_output = 1;
			break;
		case '?':
			if (optopt == 'c')
				sprintf(buffer, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				sprintf(buffer, "Unknown option `-%c'.\n", optopt);
			else
				sprintf(buffer, "Unknown option character `\\x%x'.\n", optopt);

			error(buffer);
			
		default:
			error("Unrecognized flag.");
		}
	}

	if(!tflag)
		error("Please specify the transcript file using the -t flag");

	if(!mflag)
		error("Please specify the miRNA file using the -m flag");
	
	if(!aflag)
		error("Please specify the align file using the -a flag");
		
	if(args->normalization && args->accessibility)
		error("Sorry, normalization is not available when calculating accessibility scores");
		
	args->start_time = ctime(&mytime);
	
	snprintf(args->param_info,	sizeof args->param_info,
							"# Complementarity threshold: %.1f\n"
							"# Accessibility threshold: %.1f\n"
							"# Seed region from nucleotides %u to %u\n"
							"# Central region from nucleotides %u to %u\n"
							"# Scoring schema - Match (%.1f) Mismatch(%.1f) Gaps(%.1f) Wobbles(%.1f)\n" 
							"# Accessibility: %s\tNormalization: %s\tHuman readable output: %s\n\n", 
							args->c_threshold, args->a_threshold, 
							s->seed_start, s->seed_stop, 
							s->cr_start, s->cr_stop,
							s->match, s->mismatch, s->gap, s->wobble,
							(args->accessibility ? "Yes" : "No"),
							(args->normalization ? "Yes" : "No"), (args->human_output ? "Yes" : "No"));
	
	args->temp_file = (char**)safe_malloc(sizeof(char*) * args->num_processors);
	for(i = 0; i < args->num_processors; i++)
		args->temp_file[i] = tempnam("output", "pine");
		
	return args;
}

void target_prediction(int current, fasta_align **aligns, dataset_t *tds, dataset_t *mds, pinetree_args* args, FILE* output_file){
	
	char *align1 = aligns[current]->target_seq;
	char *align2 = aligns[current]->miRNA_seq;
	int i = aligns[current]->target_id;
	int j = aligns[current]->miRNA_id;
	int start = aligns[current]->start;
 	int length = strlen(mds->sequences[j]) - count_occurrences(align1, '-');
	
	float cscore = align_score(align1, align2, args->sschema);
	
	if(cscore <= args->c_threshold){
		double ascore = calculate_accessibility("RNAup", tds->sequences[i], start, length);
		
		/*FIXME*/
		if(ascore == INFINITY)
			ascore = 25;
		
		if(ascore <= args->a_threshold){
		
			char* align = alignment_string(align1, align2, args->sschema);
			char* reg_mechanism = mechanism(align1, align2, args->sschema);
			
			if(args->human_output){
				fprintf(output_file, "target id: %s\n", tds->ids[i]);
				if(tds->annotations)
					fprintf(output_file, "target info: %s\n", tds->annotations[i] ? tds->annotations[i] : "N/A");
					
				fprintf(output_file, "miRNA id: %s\n", mds->ids[j]);
				fprintf(output_file, "Reg. mechanism: %s\n", reg_mechanism);
				fprintf(output_file, "Complementarity score: %.1f\n", cscore);
				fprintf(output_file, "Accessibility score: %.1f\n", ascore);
				fprintf(output_file, "miRNA:\t%s\n", align2);
				fprintf(output_file, "\t%s\n", align);
				fprintf(output_file, "target:\t%s\n", align1);
				fprintf(output_file, "#\n");
			}	
			else {
				fprintf(output_file, "%s,%s,%.1f,%.1f,%s,%s", tds->ids[i], mds->ids[j], cscore, ascore, align, reg_mechanism);
			
				if(tds->annotations)
					fprintf(output_file, ",%s", tds->annotations[i] ? tds->annotations[i] : "N/A");
				fprintf(output_file, "\n");
			}
			
			safe_free(align);
		}
	}
}

void simple_target_prediction(fasta_info *info, pinetree_args* args, dataset_t *tds, dataset_t *mds){
	
	uint current;
	fasta_align **aligns = info->aligns;
	
	#pragma omp parallel num_threads(args->num_processors)
	{
		/*FIXME: use O_EXCL */
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");
		
		#pragma omp for 
		for(current = 0; current < info->count; current++){
			char *align1 = aligns[current]->target_seq;
			char *align2 = aligns[current]->miRNA_seq;
			int i = aligns[current]->target_id;
			int j = aligns[current]->miRNA_id;
			float cscore = align_score(align1, align2, args->sschema);
			
			if(cscore <= args->c_threshold){
				char* align = alignment_string(align1, align2, args->sschema);
				char* reg_mechanism = mechanism(align1, align2, args->sschema);
				
				if(args->human_output){
					fprintf(output_file, "target id: %s\n", tds->ids[i]);
					if(tds->annotations)
						fprintf(output_file, "target info: %s\n", tds->annotations[i] ? tds->annotations[i] : "N/A");
						
					fprintf(output_file, "miRNA id: %s\n", mds->ids[j]);
					fprintf(output_file, "Reg. mechanism: %s\n", reg_mechanism);
					fprintf(output_file, "Complementarity score: %.1f\n", cscore);
					fprintf(output_file, "miRNA:\t%s\n", align2);
					fprintf(output_file, "\t%s\n", align);
					fprintf(output_file, "target:\t%s\n", align1);
					fprintf(output_file, "#\n");
				}	
				else {
					fprintf(output_file, "%s,%s,%.1f,%s,%s", tds->ids[i], mds->ids[j], cscore, align, reg_mechanism);
				
					if(tds->annotations)
						fprintf(output_file, ",%s", tds->annotations[i] ? tds->annotations[i] : "N/A");
					fprintf(output_file, "\n");
				}
				
				safe_free(align);
			}
		}
		
		safe_fclose(output_file);
	}
}

void normalized_prediction(fasta_info *info, pinetree_args* args, dataset_t *tds, dataset_t *mds){
	
	uint i, current;
	fasta_align **aligns = info->aligns;
	float cdiff = 0.0;
	float cmin = MAX_CSCORE, cmax = MIN_CSCORE;
	
	#pragma omp parallel for reduction(max:cmax) reduction(min:cmin) num_threads(args->num_processors) 
	for(current = 0; current < info->count; current++){
		char *align1 = aligns[current]->target_seq;
		char *align2 = aligns[current]->miRNA_seq;
		float cscore = align_score(align1, align2, args->sschema);
		aligns[current]->cscore = cscore;
		
		if(cscore < cmin)
			cmin = cscore;
			
		if(cscore > cmax)
			cmax = cscore;
	}
	
	cdiff = cmax-cmin;
	
	#pragma omp parallel num_threads(args->num_processors)
	{
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");
		
		#pragma omp for 
		for(i = 0; i < info->count; i++){
			char *align1 = aligns[i]->target_seq;
			char *align2 = aligns[i]->miRNA_seq;
			int target_id = aligns[i]->target_id;
			int miRNA_id = aligns[i]->miRNA_id;
			float cnorm = (aligns[i]->cscore - cmin)/cdiff;
			
			if(cnorm <= args->nc_threshold){
				char* align = alignment_string(align1, align2, args->sschema);
				char* reg_mechanism = mechanism(align1, align2, args->sschema);
				fprintf(output_file, "%s,%s,%.5f,%s,%s", tds->ids[target_id], mds->ids[miRNA_id], cnorm, align, reg_mechanism);
				
				if(tds->annotations)
					fprintf(output_file, ",%s", tds->annotations[target_id] ? tds->annotations[target_id] : "N/A");
				fprintf(output_file, "\n");
				
				safe_free(align);
			}
		}
		
		safe_fclose(output_file);
	}
}

int main(int argc, char **argv){
	uint j;
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	
	fasta_info *info = process_alignment(args->align_file);
	char *header = "gene,miRNA,cscore,ascore,align,reg_mechanism\n";
	
	if(args->annotation_file){
		header = "gene,miRNA,cscore,ascore,align,reg_mechanism,annotation\n";
		annotate_targets(tds, args->annotation_file);
	}
	
	/*We have to use this fork method instead of using OpenMP beacuse RNAup doesn't 
	 * like when we use parallellize code in the outer execution... */
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
				FILE *output_file = safe_fopen(args->temp_file[id], "w");
				
				for (i = lower; i < upper; i++) 
					target_prediction(i,info->aligns,tds,mds,args,output_file);
				
				safe_fclose(output_file);
				exit(0);
			}
		}
		while ((wpid = wait(&status)) > 0); /* Wait for children to exit. */
	}
	else {
		if(args->normalization)
			normalized_prediction(info, args, tds, mds);
		else /* default mode */
			simple_target_prediction(info, args, tds, mds);
	}
	
	file_joiner(args, header);

	clean_alignments(info);
	destroy_dataset(tds);
	destroy_dataset(mds);
	destroy_score_schema(args->sschema);
	
	for(j = 0; j < args->num_processors; j++)
		safe_free(args->temp_file[j]);
		
	safe_free(args->temp_file);
	safe_free(args);

	return 0;
}



