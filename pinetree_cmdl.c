#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "pinetree_cmdl.h"

pinetree_args *initialize_args(){
	pinetree_args *args = (pinetree_args*)safe_malloc(sizeof(pinetree_args));
	
	args->sschema = create_score_schema(MATCH, MISMATCH, GAP, GU, SEED_START, SEED_STOP, CR_START, CR_STOP);
	args->start_time = get_system_time();
	args->num_processors = NUM_PROCESSORS;
	args->transcript_file = NULL;
	args->mirna_file = NULL;
	args->output_file = "pinetree";
	args->annotation_file = NULL;
	args->c_threshold = NC_THRESHOLD;
	args->a_threshold = NA_THRESHOLD;
	args->evalue = EVALUE;
	args->accessibility = 0;
	args->normalization = 1;
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
		else if(STRMATCH("SEED_REGION_START", parameter))
			sscanf(value, "%u", &(s->seed_start));
		else if(STRMATCH("SEED_REGION_STOP", parameter))
			sscanf(value, "%u", &(s->seed_stop));
		else if(STRMATCH("CENTRAL_REGION_START", parameter))
			sscanf(value, "%u", &(s->cr_start));
		else if(STRMATCH("CENTRAL_REGION_STOP", parameter))
			sscanf(value, "%u", &(s->cr_stop));
		else if(STRMATCH("EVALUE", parameter))
			sscanf(value, "%u", &(args->evalue));
		else if(STRMATCH("COMP_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->c_threshold));
		else if(STRMATCH("ACC_THRESHOLD", parameter))
			sscanf(value, "%f", &(args->a_threshold));
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
		BOLDWHITE "-b, --accessibility\n" RESET 
		"\tturns on accessibility mode\n"
		BOLDWHITE "-h, --help\n" RESET 
		"\tdetailed information about the command-line options\n"
		BOLDWHITE "-n, --no_normalization\n" RESET 
		"\tturns off normalization mode \n"
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
		BOLDWHITE "-a, --a_cut [value]\n" RESET 
		"\tcutoff value for the accessibility criterion\n"
		BOLDWHITE "-A, --annotation [annotation file]\n" RESET 
		"\tprovides information about the transcripts functions\n"
		BOLDWHITE "-b, --accessibility\n" RESET 
		"\tturns on accessibility mode\n"
		BOLDWHITE "-c, --c_cut [value]\n" RESET 
		"\tcutoff value for the complementarity criterion\n"
		BOLDWHITE "-C, --config [file]\n" RESET 
		"\tindicates where the path of the CONFIG file\n"
		BOLDWHITE "-e, --evalue [value]\n" RESET 
		"\tcutoff value for the e-value parameter of FASTA\n"
		BOLDWHITE "-h, --help\n" RESET 
		"\tdetailed information about the command-line options\n"
		BOLDWHITE "-n, --no_normalization\n" RESET 
		"\tturns off normalization mode \n"
		BOLDWHITE "-o, --output [file name]\n" RESET 
		"\tname of the file where the output os stored (e.g -o pinetree)\n"
		BOLDWHITE "-p --processors\n" RESET 
		"\tnumber of processors\n"
		BOLDWHITE "-P, --pretty\n" RESET 
		"\tenable human readable output\n"
		BOLDWHITE "-v, --version\n" RESET 
		"\tprints the current version of the program\n"
		"\n"
		BOLDWHITE "Scoring schema parameters:\n" RESET 
		BOLDWHITE "-g, --gap [value]\n" RESET 
		"\tvalue for the gap penalty\n"
		BOLDWHITE "-u, --match [value]\n" RESET 
		"\tvalue for the match\n"
		BOLDWHITE "-w, --wobble [value]\n" RESET 
		"\tvalue for the G:U wobble\n"
		BOLDWHITE "-y, --mismatch [value]\n" RESET 
		"\tvalue for the mismatch\n"
		BOLDWHITE "-x, --central_start [value]\n" RESET 
		"\tindicates where the central region starts\n"
		BOLDWHITE "-X, --central_stop [value]\n" RESET 
		"\tindicates where the central region ends\n"
		BOLDWHITE "-z, --seed_start [value]\n" RESET 
		"\tindicates where the seed region starts\n"
		BOLDWHITE "-Z, --seed_stop [value]\n" RESET 
		"\tindicates where the seed region ends\n"
		);
	exit(0);
}

pinetree_args* read_cml_arguments(int argc, char **argv){
	int i, c;
	boolean tflag = 0, mflag = 0;
	boolean aflag = 0, cflag = 0;
	char *config_file = NULL;
	pinetree_args *args = initialize_args();
	score_t *s = args->sschema;

	while (1){
      
		int option_index = 0;
		static struct option long_options[] = {
			{"annotation", required_argument, 0, 'A'},
			{"a_cut", required_argument, 0, 'a'},
			{"accessibility", no_argument, 0, 'b'},
			{"central_start", required_argument, 0, 'x'},
			{"central_stop", required_argument, 0, 'X'},
			{"config", required_argument, 0, 'C'}, 
			{"c_cut", required_argument, 0, 'c'},
			{"cnorm_cut", required_argument, 0, 'd'},
			{"evalue", required_argument, 0, 'e'},
			{"gap",  required_argument, 0, 'g'},
			{"help", no_argument, 0, 'h'},
			{"match",  required_argument, 0, 'u'},
			{"mismatch",  required_argument, 0, 'y'},
			{"mirna",  required_argument, 0, 'm'},
			{"no_normalization", no_argument, 0, 'n'},
			{"output",  required_argument, 0, 'o'},
			{"pretty",  no_argument, 0, 'P'},
			{"processors",  required_argument, 0, 'p'},
			{"seed_start", required_argument, 0, 'z'},
			{"seed_stop", required_argument, 0, 'Z'},
			{"target",  required_argument, 0, 't'},
			{"version",  no_argument, 0, 'v'},
			{"wobble",  required_argument, 0, 'w'},
			{0, 0, 0, 0}
		};

		c = getopt_long (argc, argv, "a:A:bc:C:e:f:g:hm:no:p:Pt:u:vw:x:X:y:z:Z", 
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;
      
		switch (c){
		case 'a':
			args->a_threshold = atof(optarg);
			aflag = 1;
			break;
		case 'A':
			args->annotation_file = optarg;
			break;
		case 'b':
			args->accessibility = 1;
			break;
		case 'c':
			args->c_threshold = atof(optarg);
			cflag = 1;
			break;
		case 'C':
			config_file = optarg;
			initialize_parameters(config_file, args);
			break;
		case 'e':
			args->evalue = atoi(optarg);
			break;
		case 'g':
			s->gap = atof(optarg);
			break;
		case 'h':
			print_usage();
		case 'm':
			mflag = 1;
			args->mirna_file = optarg;
			break;
		case 'n':
			args->normalization = 0;
			if(!aflag)
				args->a_threshold = A_THRESHOLD;
			if(!cflag)
				args->c_threshold = C_THRESHOLD;
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
		case 'u':
			s->match = atof(optarg);
			break;
		case 'v':
			print_version();
			exit(0);
		case 'w':
			s->wobble = atof(optarg);
			break;
		case 'x':
			s->cr_start = atoi(optarg);
			break;
		case 'X':
			s->cr_stop = atoi(optarg);
			break;
		case 'y':
			s->mismatch = atof(optarg);
			break;
		case 'z':
			s->seed_start = atoi(optarg);
			break;
		case 'Z':
			s->seed_stop = atoi(optarg);
			break;
		case '?':
			exit(-1);
		default:
			print_usage();
		}
	}

	if(!tflag || !mflag)
		common_usage();
		
	snprintf(args->header, sizeof args->header,
							"gene,miRNA,cscore%s,align,reg_mechanism%s\n",
							(args->accessibility ? ",ascore" : ""),
							(args->annotation_file ? ",annotation" : ""));
	
	snprintf(args->param_info,	sizeof args->param_info,
							"# e-value threshold: %d\n"
							"# Complementarity threshold: %.*f\n"
							"# Accessibility threshold: %.*f\n"
							"# Seed region from nucleotides %u to %u\n"
							"# Central region from nucleotides %u to %u\n"
							"# Scoring schema - Match(%.1f) Mismatch(%.1f) Gaps(%.1f) Wobbles(%.1f)\n" 
							"# Accessibility: %s\tNormalization: %s\tHuman readable output: %s\n\n", 
							args->evalue,
							args->normalization ? 5 : 1, args->c_threshold, 
							args->normalization ? 5 : 1, args->a_threshold, 
							s->seed_start, s->seed_stop, 
							s->cr_start, s->cr_stop,
							s->match, s->mismatch, s->gap, s->wobble,
							(args->accessibility ? "Yes" : "No"),
							(args->normalization ? "Yes" : "No"), (args->human_output ? "Yes" : "No"));
	
	args->temp_file = (char**)safe_malloc(sizeof(char*) * args->num_processors);
	for(i = 0; i < args->num_processors; i++){
		char *temp_file = create_unique_file("pinetree_XXXXXX");
		args->temp_file[i] = temp_file;
	}
		
	return args;
}
