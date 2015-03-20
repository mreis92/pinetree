#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <getopt.h>
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
	args->output_file = "pinetree";
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
	
	safe_fclose(config_file);
}


void common_usage(){
	print_version();
	
	fprintf(stderr,
		BOLDWHITE "USAGE:\n" RESET 
		"\t bash pinetree.sh -t [transcript file] -m [mirna file] [-options]\n\n"
		
		BOLDWHITE "Common options:\n" RESET 
		BOLDWHITE "-A, --annotation [annotation file]\n" RESET 
		"\tprovides information about the transcripts functions\n"
		BOLDWHITE "-h, --help\n" RESET 
		"\tdetailed information about the command-line options\n"
		BOLDWHITE "-o, --output [file name]\n" RESET 
		"\tname of the file where the output os stored (e.g -o pinetree)\n"
		BOLDWHITE "-p --processors\n" RESET 
		"\tnumber of processors\n"
		BOLDWHITE "-P, --pretty\n" RESET 
		"\tenable human readable output\n"
		);
	exit(0);
}

void print_usage(){
	print_version();
	
	fprintf(stderr,
		BOLDWHITE "USAGE:\n" RESET 
		"\t bash pinetree.sh -t [transcript file] -m [mirna file] [-options]\n"
		"\n"
		BOLDWHITE "Available options:\n" RESET 
		BOLDWHITE "-A, --annotation [annotation file]\n" RESET 
		"\tprovides information about the transcripts functions\n"
		BOLDWHITE "-c, --c_cut [value]\n" RESET 
		"\tcutoff value for the complementarity criterion\n"
		BOLDWHITE "-C, --config [file]\n" RESET 
		"\tindicates where the path of the CONFIG file\n"
		BOLDWHITE "-e, --e_cut [value]\n" RESET 
		"\tcutoff value for the anti-target criterion\n"
		BOLDWHITE "-h, --help\n" RESET 
		"\tdetailed information about the command-line options\n"
		BOLDWHITE "-o, --output [file name]\n" RESET 
		"\tname of the file where the output os stored (e.g -o pinetree)\n"
		BOLDWHITE "-p --processors\n" RESET 
		"\tnumber of processors\n"
		BOLDWHITE "-P, --pretty\n" RESET 
		"\tenable human readable output\n"
		BOLDWHITE "-v, --version\n" RESET 
		"\tprints the current version of the program\n"
		"\n"
		BOLDWHITE "Anti-target parameters:\n" RESET 
		BOLDWHITE "-E, --errors [value]\n" RESET 
		"\tvalue for the allowed number of errors\n"
		BOLDWHITE "-M, --markov [value]\n" RESET 
		"\tvalue for the Markov order of the modules\n"
		);
	exit(0);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	char buffer[BUFSIZE];
	int i, c;
	int tflag = 0, mflag = 0;
	char *config_file = NULL;
	pinetree_args *args = initialize_args();
	
	while (1){
      
		int option_index = 0;
		static struct option long_options[] = {
			{"annotation", required_argument, 0, 'A'},
			{"config", required_argument, 0, 'C'}, 
			{"e_cut", required_argument, 0, 'e'},
			{"errors", required_argument, 0, 'E'},
			{"help", no_argument, 0, 'h'},
			{"mirna",  required_argument, 0, 'm'},
			{"markov", required_argument, 0, 'M'},
			{"output",  required_argument, 0, 'o'},
			{"pretty",  no_argument, 0, 'P'},
			{"processors",  required_argument, 0, 'p'},
			{"target",  required_argument, 0, 't'},
			{0, 0, 0, 0}
		};

		c = getopt_long (argc, argv, "A:C:e:E:hm:M:o:p:Pt:", 
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;
      
		switch (c){
		case 'A':
			args->annotation_file = optarg;
			break;
		case 'C':
			config_file = optarg;
			initialize_parameters(config_file, args);
			break;
		case 'e':
			args->e_threshold = atof(optarg);
			break;
		case 'E':
			args->num_errors = atoi(optarg);
			break;
		case 'h':
			print_usage();
			exit(0);
		case 'm':
			mflag = 1;
			args->mirna_file = optarg;
			break;
		case 'M':
			args->markov_order = atoi(optarg);
			break;
		case 'o':
			args->output_file = optarg;
			break;
		case 'p':
			args->num_processors = atoi(optarg);
			break;
		case 'P':
			args->human_output = 1;
			break;
		case 't':
			tflag = 1;
			args->transcript_file = optarg;
			break;
		default:
			print_usage();
			exit(0);
		}
	}
	
	if(!tflag || !mflag)
		common_usage();
		
	snprintf(args->header, sizeof args->header,
							"gene,miRNA,escore%s\n",
							(args->annotation_file ? ",annotation" : ""));
		
	snprintf(args->param_info,	sizeof args->param_info,
							"# Evolutionary threshold: %.1f\n"
							"# Markov order: %u\n"
							"# Number of errors allowed: %u\n"
							"# Human readable output: %s\n\n",
							args->e_threshold, 
							args->markov_order, 
							args->num_errors,
							(args->human_output ? "Yes" : "No"));
	
	args->temp_file = (char**)safe_malloc(sizeof(char*) * args->num_processors);
	for(i = 0; i < args->num_processors; i++){
		char *temp_file = create_unique_file("/tmp/pinetree_XXXXXX");
		args->temp_file[i] = temp_file;
	}
	
	return args;
}

int main(int argc, char **argv){
	uint proc;
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	evo_info_t* evo_info = generate_models(tds, tds, mds, args->markov_order, args->num_processors);
	
	if(args->annotation_file)
		annotate_targets(tds, args->annotation_file);
	
	
	#pragma omp parallel
	{
		uint i, j;
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");

		#pragma omp for 
		for (i = 0; i < tds->seqn; i++) {
			for(j = 0; j < mds->seqn; j++){
				ldouble escore = calculate_escore(evo_info, j, i);
				
				if(args->human_output){
					fprintf(output_file, "target id: %s\n", tds->ids[i]);
					if(tds->annotations)
						fprintf(output_file, "target info: %s\n", tds->annotations[i] ? tds->annotations[i] : "N/A");
						
					fprintf(output_file, "miRNA id: %s\n", mds->ids[j]);
					fprintf(output_file, "Evolutionary score: %Lg\n", escore);
					fprintf(output_file, "#\n");
				}	
				else {
					fprintf(output_file, "%s,%s,%Lg", tds->ids[i], mds->ids[j], escore);
					
					if(tds->annotations)
						fprintf(output_file, ",%s", tds->annotations[i] ? tds->annotations[i] : "N/A");
					fprintf(output_file, "\n");
				}
			}
		}
		
		safe_fclose(output_file);
	}
	
	file_joiner(args);
	
	cleanup(evo_info); 
	destroy_dataset(tds);
	destroy_dataset(mds);
	
	for(proc = 0; proc < args->num_processors; proc++)
		safe_free(args->temp_file[proc]);
		
	safe_free(args->temp_file);
	safe_free(args);

	return 0;
}



