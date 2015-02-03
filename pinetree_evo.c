#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <sys/wait.h>
#include <omp.h>
#include <unistd.h>

#include "annotation.h"
#include "constants.h"
#include "dataset.h"
#include "evo.h"
#include "fasta.h"
#include "pinetree_utils.h"
#include "util.h"

pinetree_args *initialize_args(){
	pinetree_args *args = (pinetree_args*) safe_malloc(sizeof(pinetree_args));
	args->start_time = get_system_time();
	args->num_processors = NUM_PROCESSORS;
	args->transcript_file = NULL;
	args->mirna_file = NULL;
	args->output_file = NULL;
	args->annotation_file = NULL;
	args->e_threshold = E_THRESHOLD;
	args->num_errors = NUM_ERRORS;
	args->markov_order = MARKOV_ORDER;
	args->human_output = 0;
	
	return args;
}

void initialize_parameters(char* filename, pinetree_args* args){
	char buffer[BUFSIZE];
	char *parameter, *value;
	uint i;

	FILE *config_file = safe_fopen(filename, "r");

	while(fgets (buffer, BUFSIZE, config_file) != NULL){
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
	
	snprintf(args->param_info,	sizeof args->param_info,
							"# Evolutionary threshold: %.1f\n"
							"# Markov order: %u\n"
							"# Number of errors allowed: %u\n",
							args->e_threshold, 
							args->markov_order, 
							args->num_errors);
	
	args->temp_file = (char**)safe_malloc(sizeof(char*) * args->num_processors);
	for(i = 0; i < args->num_processors; i++)
		args->temp_file[i] = tempnam("output", "pine");
	
	safe_fclose(config_file);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	char buffer[BUFSIZE];
	int c;
	int tflag = 0, mflag = 0;
	char *config_file = NULL;
	pinetree_args *args = initialize_args();
	
	opterr = 0;
	
	while ((c = getopt (argc, argv, "C:t:m:o:n:A:")) != -1){
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
	
	return args;
}

int main(int argc, char **argv){
	uint proc;
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	evo_info_t* evo_info = generate_models(tds, tds, mds, args->markov_order, args->num_processors);
	
	char* header = "gene,miRNA,escore\n";
	
	if(args->annotation_file){
		header = "gene,miRNA,escore,annotation\n";
		annotate_targets(tds, args->annotation_file);
	}
	
	#pragma omp parallel
	{
		uint i, j;
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");

		#pragma omp for 
		for (i = 0; i < tds->seqn; i++) {
			for(j = 0; j < mds->seqn; j++){
				long double escore = calculate_escore(evo_info, j, i);
				fprintf(output_file, "%s,%s,%Lg", tds->ids[i], mds->ids[j], escore);
				
				if(tds->annotations)
					fprintf(output_file, ",%s", tds->annotations[i] ? tds->annotations[i] : "N/A");
				fprintf(output_file, "\n");
			}
		}
		
		safe_fclose(output_file);
	}
	
	file_joiner(args, header);
	
	cleanup(evo_info); 
	destroy_dataset(tds);
	destroy_dataset(mds);
	
	for(proc = 0; proc < args->num_processors; proc++)
		safe_free(args->temp_file[proc]);
		
	safe_free(args->temp_file);
	safe_free(args);

	return 0;
}



