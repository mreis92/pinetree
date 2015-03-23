#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "annotation.h"


/* Annotates an existing dataset with transcript specific information */
void annotate_targets(dataset_t *tds, StrMap *id_map, char* filename){
	char gene_id[BUFSIZE], description[LONGBUF];
	char id[IDSIZE];
	FILE *annotation_file = safe_fopen(filename, "r");
	
	tds->annotations = (char**)safe_calloc(tds->seqn, sizeof(char*));
	
	while (fscanf(annotation_file,"%[^\t]\t%[^\n]\n", gene_id, description) != EOF){
		sm_get(id_map, gene_id, id, IDSIZE);
		tds->annotations[atoi(id)] = strdup(description);
	}
		
	safe_fclose(annotation_file);
	sm_delete(id_map);
}
