#include <stdio.h>
#include <string.h>

#include "accessibility.h"

double calculate_accessibility(char* mode, char *gene_seq, int start, int align_length){
	char *seq;
	int init_pos, end_pos, up_flank, down_flank, region_length;
	int gene_len = strlen(gene_seq);

	if(align_length+start>gene_len)
		align_length = gene_len-start+1;
	
	/* Testing the limits for flanks and up/downstream parameters */
	if (start-FLANK > 0){
		init_pos =  (start-1)-FLANK;
		up_flank = FLANK;
	}
	else {
		init_pos = 0;
		up_flank = start-1;
	}

	if ((start+align_length+FLANK) < gene_len) {
		end_pos = (start-1)+align_length+FLANK ;
		down_flank = DOWNSTREAM;
	}
	else {
		end_pos = gene_len;
		down_flank = 0; /*FIXME: there should be some calculations here */
	}
	
	/* Run acc module with this parameters */
	region_length = end_pos-init_pos+1; /* include \0 terminator */
	seq = (char*) safe_malloc(sizeof(char) * region_length);

	snprintf(seq, region_length, "%s", gene_seq + init_pos);

	if(!strcmp(mode, "RNAup"))
		return accessibility_RNAup(seq, align_length, up_flank+align_length);
	else 
		error("Unrecognized accessibility mode.");
}

double accessibility_RNAup(char *gene, int region, int end_pos){
	char 	program[] = "RNAup";
	/*char  	mode_flag[] = "--interaction_pairwise";*/
	char  	region_flag[] = "-u";
	char  	window_flag[] = "-w";
	char    threeprimereg_flag[] = "-3";
	char    fiveprimereg_flag[] = "-5";
	char  	region_length[3]; 
	char 	window_length[4];
	char    threeprimereg_length[3];
	char    fiveprimereg_length[3];
	sprintf(fiveprimereg_length, "%d", UPSTREAM);
	sprintf(threeprimereg_length, "%d", DOWNSTREAM);
	sprintf(region_length, "%d", region+UPSTREAM+DOWNSTREAM);
	sprintf(window_length, "%d", region);
	char* 	argv[] = { &program[0], &region_flag[0], &region_length[0], &window_flag[0], &window_length[0], &fiveprimereg_flag[0],  &fiveprimereg_length[0], &threeprimereg_flag[0], &threeprimereg_length[0],  NULL };
   	int 	argc    = (int)(sizeof(argv) / sizeof(argv[0])) - 1;

	return RNAup(argc, argv, gene, region, end_pos);
}

