#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/wait.h>
#include <unistd.h>
#include "accessibility.h"
#include "align.h"
#include "annotation.h"
#include "dataset.h"
#include "constants.h"
#include "fasta.h"
#include "util.h"

typedef struct pinetree_args {
	score_t *sschema;
	uint num_processors;
	char* transcript_file;
	char* mirna_file; 
	char* align_file; 
	char* annotation_file;
	char* output_file;
	float c_threshold;
	float a_threshold;
} pinetree_args;

typedef struct fasta_align{
	char* miRNA_seq;
	char* target_seq;
	int target_id;
	int miRNA_id;
	int start;
} fasta_align;

typedef struct fasta_info{
	fasta_align **aligns;
	int count;
} fasta_info;

pinetree_args *initialize_args(){
	pinetree_args *args = (pinetree_args*) malloc(sizeof(pinetree_args));
	
	args->sschema = create_score_schema(MATCH, MISMATCH, GAP, GU, SEED_PENALTY, SEED_START, SEED_STOP, CR_START, CR_STOP);
	args->num_processors = NUM_PROCESSORS;
	args->transcript_file = NULL;
	args->mirna_file = NULL;
	args->align_file = NULL;
	args->output_file = NULL;
	args->annotation_file = NULL;
	args->c_threshold = C_THRESHOLD;
	args->a_threshold = A_THRESHOLD;
	
	return args;
}

void file_joiner(pinetree_args* args, char* header){
	uint i;
	FILE *output_file;
	char buffer[BUFSIZE] = "";
	
	snprintf(buffer, BUFSIZE, "%s.csv", args->output_file);
	output_file = safe_fopen(buffer, "w");
	fputs(header, output_file);
  
	for (i = 0; i < args->num_processors; ++i) {
		FILE *file;
		char line[BUFSIZE*2]; 

		snprintf(buffer, BUFSIZE, "%s_tmp_%d.csv", args->output_file, i);
		file = safe_fopen(buffer, "r");

		while(fgets(line, BUFSIZE*2, file) != NULL){
			fputs(line,output_file); 
		}

		safe_fclose(file);
		safe_remove(buffer);
	}
	
	safe_fclose(output_file);
}

void initialize_parameters(char* filename, pinetree_args* args){
	char buffer[50];
	char *parameter, *value;

	FILE *config_file = safe_fopen(filename, "r");

	while(fgets (buffer, 50, config_file) != NULL){
		parameter = strtok(buffer, "=");
		value = strtok(NULL, "");

		if(STRMATCH("NUM_PROC", parameter))
			sscanf(value, "%u", &(args->num_processors));
		else if(STRMATCH("GAP", parameter))
			sscanf(value, "%f", &(args->sschema->gap));
		else if(STRMATCH("MATCH", parameter))
			sscanf(value, "%f", &(args->sschema->match));
		else if(STRMATCH("MISMATCH", parameter))
			sscanf(value, "%f", &(args->sschema->mismatch));
		else if(STRMATCH("GU_WOBBLE", parameter))
			sscanf(value, "%f", &(args->sschema->wobble));
		else if(STRMATCH("SEED_PENALTY", parameter))
			sscanf(value, "%f", &(args->sschema->seed_penalty));
		else if(STRMATCH("SEED_REGION_LOWER_BOUND", parameter))
			sscanf(value, "%u", &(args->sschema->seed_start));
		else if(STRMATCH("SEED_REGION_UPPER_BOUND", parameter))
			sscanf(value, "%u", &(args->sschema->seed_stop));
		else if(STRMATCH("CENTRAL_REGION_LOWER_BOUND", parameter))
			sscanf(value, "%u", &(args->sschema->cr_start));
		else if(STRMATCH("CENTRAL_REGION_UPPER_BOUND", parameter))
			sscanf(value, "%u", &(args->sschema->cr_stop));
		else if(STRMATCH("COMP_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->c_threshold));
		else if(STRMATCH("ACC_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->a_threshold));
	}
	
	safe_fclose(config_file);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	char buffer[100];
	int c;
	int tflag = 0, mflag = 0, aflag = 0;
	char *config_file = NULL;
	
	pinetree_args *args = initialize_args();

	opterr = 0;
	
	while ((c = getopt (argc, argv, "C:t:m:a:o:n:A:")) != -1){
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

	return args;
}

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

fasta_info *process_alignment(char* filename){
	fasta_info *info;
	fasta_align **aligns;
	int num_entries;
	int miRNA_id, target_id;
	char miRNA_seq[BUFSIZE], target_seq[BUFSIZE];
	int start;
	
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

void target_prediction(int current, fasta_align **aligns, dataset_t *tds, dataset_t *mds, pinetree_args* args, FILE* output_file){
	
	char *align1 = aligns[current]->target_seq;
	char *align2 = aligns[current]->miRNA_seq;
	int i = aligns[current]->target_id;
	int j = aligns[current]->miRNA_id;
	int start = aligns[current]->start;
 	int length = strlen(mds->sequences[j]) - count_occurrences(align1, '-');
	
	float gscore = align_score(align1, align2, args->sschema);
	
	if(gscore <= args->c_threshold){
		double ascore = calculate_accessibility("RNAup", tds->sequences[i], start, length);
		
		if(ascore <= args->a_threshold){
		
			char* align = alignment_string(align1, align2, args->sschema);
			char* reg_mechanism = mechanism(align1, align2, args->sschema);
			fprintf(output_file, "%s,%s,%.1f,%.1f,%s,%s", tds->ids[i], mds->ids[j], gscore, ascore, align, reg_mechanism);
			if(tds->annotations)
				fprintf(output_file, ",%s", tds->annotations[i] ? tds->annotations[i] : "N/A");
			fprintf(output_file, "\n");
			
			safe_free(align);
		}
	}
}

int main(int argc, char **argv){
	uint id;
	int status;
	pid_t wpid;
	
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	
	fasta_info *info = process_alignment(args->align_file);
	char *header = "gene,miRNA,gscore,ascore,align,reg_mechanism\n";
	
	if(args->annotation_file){
		header = "gene,miRNA,gscore,ascore,align,reg_mechanism,annotation\n";
		annotate_targets(tds, args->annotation_file);
	}
	
	int batch = info->count/args->num_processors;
	int remainder = info->count%args->num_processors;
	
	/*uint i;
	FILE *output_file;
	uint lower = id*batch;
	uint upper = (id==args->num_processors-1) ? (id+1)*batch + remainder : (id+1)*batch;
	char buffer[BUFSIZE] = "";
	snprintf(buffer, BUFSIZE, "%s_tmp_%d.csv", args->output_file, 0);
	output_file = safe_fopen(buffer, "w");
			
	for(i = 0; i < info->count; i++){
		target_prediction(i,info->aligns,tds,mds,args,output_file);
	}*/
	
	/* Start children. */
	for (id = 0; id < args->num_processors; ++id) {
		if ((wpid = fork()) < 0) 
			printf("Error in fork\n");
		else if (wpid == 0) {
			uint i;
			FILE *output_file;
			uint lower = id*batch;
			uint upper = (id==args->num_processors-1) ? (id+1)*batch + remainder : (id+1)*batch;
			char buffer[BUFSIZE] = "";
			
			snprintf(buffer, BUFSIZE, "%s_tmp_%d.csv", args->output_file, id);
			output_file = safe_fopen(buffer, "w");
		
			for (i = lower; i < upper; i++) {
				target_prediction(i,info->aligns,tds,mds,args,output_file);
			}
			
			safe_fclose(output_file);
			exit(0);
		}
	}

	/* Wait for children to exit. */
	while ((wpid = wait(&status)) > 0);
	
	file_joiner(args,header);

	clean_alignments(info);
	destroy_dataset(tds);
	destroy_dataset(mds);
	destroy_score_schema(args->sschema);
	safe_free(args);

	return 0;
}



