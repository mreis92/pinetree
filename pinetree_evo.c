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
#include "pinetree_evo_cmdl.h"
#include "pinetree_utils.h"
#include "util.h"

int main(int argc, char **argv){
	pinetree_args* args = read_cml_arguments(argc, argv);
	
	dataset_t *tds = parse_fasta(args->transcript_file);
	dataset_t *mds = parse_fasta(args->mirna_file);
	evo_info_t* evo_info = generate_models(tds, tds, mds, args->markov_order, args->num_processors);
	
	StrMap *target_map = map_ids(tds);
	
	if(args->annotation_file)
		annotate_targets(tds, target_map, args->annotation_file);
	
	sm_delete(target_map);
	
	#pragma omp parallel num_threads(args->num_processors)
	{
		uint i, j;
		FILE *output_file = safe_fopen(args->temp_file[omp_get_thread_num()], "w");

		#pragma omp for 
		for (i = 0; i < tds->seqn; i++) {
			for(j = 0; j < mds->seqn; j++){
				float affinity = calculate_affinity(evo_info, j, i);
				
				if(args->human_output){
					fprintf(output_file, "target id: %s\n", tds->ids[i]);
					if(tds->annotations)
						fprintf(output_file, "target info: %s\n", tds->annotations[i] ? tds->annotations[i] : "N/A");
						
					fprintf(output_file, "miRNA id: %s\n", mds->ids[j]);
					fprintf(output_file, "miRNA:mRNA affinity: %f\n", affinity);
					fprintf(output_file, "#\n");
				}	
				else {
					fprintf(output_file, "%s,%s,%f", tds->ids[i], mds->ids[j], affinity);
					
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
		
	safe_free(args);

	return 0;
}



