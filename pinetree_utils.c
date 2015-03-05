#include <stdio.h>
#include <time.h>

#include "pinetree_utils.h"

/* Prints the current version of the program */
void print_version(){
	fprintf(stderr, GREEN 
	"\t\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	"\t\t  PINETREE: Version %.1f\n"
	"\t\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	,RESET, VERSION);
}

/* Joins files that resulted from parallel execution in a single file */
void file_joiner(pinetree_args* args){
	uint i;
	FILE *output_file;
	char buffer[BUFSIZE] = "";

	snprintf(buffer, BUFSIZE, "%s.%s", args->output_file, (args->human_output ? "txt" : "csv"));
	output_file = safe_fopen(buffer, "w");
	
	fprintf(output_file, "# Start time of execution: %s", args->start_time);
	fprintf(output_file, "# Finish time of execution: %s", get_system_time());
	fprintf(output_file, "# Number of processors: %u\n\n", args->num_processors);
	
	fprintf(output_file, "# Transcript file used: %s\n", args->transcript_file);
	fprintf(output_file, "# miRNA file used: %s\n", args->mirna_file);
	
	if(args->annotation_file)
		fprintf(output_file, "# Annotation file used: %s\n", args->annotation_file);
		
	fputs(args->param_info, output_file);
	
	if(!args->human_output)
		fputs(args->header, output_file);
  
	for (i = 0; i < args->num_processors; ++i) {
		FILE *file;
		char line[LONGBUF]; 

		file = safe_fopen(args->temp_file[i], "r");

		while(fgets(line, LONGBUF, file) != NULL){
			fputs(line,output_file); 
		}

		safe_fclose(file);
		safe_remove(args->temp_file[i]); /* Deleting temporary file */
	}
	
	safe_fclose(output_file);
}
