#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "annotation.h"

/* Maps the string id into its position in the data structure */
StrMap *map_targets(dataset_t *tds){
	uint i;
	char id[IDSIZE];
	StrMap * map = sm_new(tds->seqn*3);
	
	for(i = 0; i < tds->seqn; i++){
		sprintf(id, "%d", i);
		sm_put(map, tds->ids[i], id);
	}
	
	return map;
}

/* Annotates an existing dataset with transcript specific information */
void annotate_targets(dataset_t *tds, char* filename){
	char gene_id[BUFSIZE], description[LONGBUF];
	char id[IDSIZE];
	StrMap *id_map = map_targets(tds);
	FILE *annotation_file = safe_fopen(filename, "r");
	
	tds->annotations = (char**)safe_calloc(tds->seqn, sizeof(char*));
	
	while (fscanf(annotation_file,"%[^\t]\t%[^\n]\n", gene_id, description) != EOF){
		sm_get(id_map, gene_id, id, IDSIZE);
		tds->annotations[atoi(id)] = strdup(description);
	}
		
	safe_fclose(annotation_file);
	sm_delete(id_map);
}
