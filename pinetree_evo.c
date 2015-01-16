#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/wait.h>
#include <omp.h>
#include <unistd.h>

#include "evo.h"
#include "dataset.h"
#include "constants.h"
#include "util.h"
#include "fasta.h"

typedef struct pinetree_args {
	uint num_processors;
	char* transcript_file;
	char* mirna_file; 
	char* output_file;
	char* subset_file;
	float e_threshold;
	uint num_errors;
	uint markov_order;
} pinetree_args;

void file_joiner(pinetree_args* args, char* header){
	uint i;

	char buffer[BUFSIZE] = "";
	snprintf(buffer, BUFSIZE, "%s.csv", args->output_file);
	FILE *output_file = safe_fopen(buffer, "w");
	fputs(header, output_file);
  
	for (i = 0; i < args->num_processors; ++i) {
		char line[BUFSIZE]; 

		snprintf(buffer, BUFSIZE, "%s_tmp_%d.csv", args->output_file, i);

		FILE *file = safe_fopen(buffer, "r" );

		fgets(line, BUFSIZE, file); /*read header */
		while(fgets(line, BUFSIZE, file) != NULL){
			fputs(line,output_file); 
		}

		fclose(file);
		remove(buffer); /*TODO: safe_remove */
	}
	fclose(output_file);
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
		else if(STRMATCH("NUM_ERRORS", parameter))
			sscanf(value, "%u", &(args->num_errors));
		else if(STRMATCH("MARKOV_ORDER", parameter))
			sscanf(value, "%u", &(args->markov_order));
		else if(STRMATCH("EVO_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->e_threshold));
	}
	
	fclose(config_file);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	char buffer[100];
	int c;
	int tflag = 0, mflag = 0;
	char *config_file = NULL;
	pinetree_args *args = (pinetree_args*) malloc(sizeof(pinetree_args));
	args->subset_file = NULL;
	args->num_processors = NUM_PROCESSORS;
	args->e_threshold = E_THRESHOLD;
	args->num_errors = NUM_ERRORS;
	args->markov_order = MARKOV_ORDER;
	
	opterr = 0;
	
	while ((c = getopt (argc, argv, "C:t:m:o:n:")) != -1){
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
		case 'o':
			args->output_file = optarg;
			break;
		case 'n':
			args->num_processors = atoi(optarg);
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
	
	return args;
}

int main(int argc, char **argv){
	
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	evo_info_t* evo_info = generate_models(tds, tds, mds, args->markov_order, args->num_processors);
	
	char* output_header = "gene,miRNA,escore\n";
	
	#pragma omp parallel
	{
	uint i, j;
	FILE *output_file;
	char buffer[BUFSIZE] = "";
	
	snprintf(buffer, BUFSIZE, "%s_tmp_%d.csv", args->output_file, omp_get_thread_num());
	output_file = safe_fopen(buffer, "w");
	fprintf(output_file, output_header);

	#pragma omp for 
	for (i = 0; i < tds->seqn; i++) {
		for(j = 0; j < mds->seqn; j++){
			long double escore = calculate_escore(evo_info, j, i);
			fprintf(output_file, "%s,%s,%Lg\n", tds->ids[i], mds->ids[j], escore);
		}
	}
	
	fclose(output_file);
	}
	
	file_joiner(args, output_header);
	
	cleanup(evo_info); 
	destroy_dataset(tds);
	destroy_dataset(mds);
	safe_free(args);

	return 0;
}



