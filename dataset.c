#include <stdlib.h>
#include <stdio.h>

#include "dataset.h"

/* Creates a dataset with multiple sequences, each one with an 
 * associated id */
dataset_t *create_dataset(uint seqn, char **sequences, char **ids){
	dataset_t *ds = (dataset_t *) safe_malloc(sizeof(dataset_t));
	ds->seqn = seqn;
	ds->ids = ids;
	ds->annotations = NULL;
	ds->sequences = sequences;

	return ds;
}

/* Frees the allocated memory of the data structure that corresponds 
 * to a dataset */
void destroy_dataset(dataset_t * ds){
	uint i = 0;

	for (i = 0; i < ds->seqn; i++) {
		safe_free(ds->sequences[i]);
		safe_free(ds->ids[i]);
	}
	
	if(ds->annotations){
		for (i = 0; i < ds->seqn; i++) {
			safe_free(ds->annotations[i]);
		}
		safe_free(ds->annotations);
	}

	safe_free(ds->sequences);
	safe_free(ds->ids);
	safe_free(ds);
}
