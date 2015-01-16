#ifndef DATASET_H
#define DATASET_H

#include "types.h"
#include "util.h"

typedef struct dataset_str {
	uint seqn;
	char **sequences;
	char **ids;
} dataset_t;

dataset_t *create_dataset(uint seqn, char **sequences, char **ids);
void destroy_dataset(dataset_t * ds);

#endif
