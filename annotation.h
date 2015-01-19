#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <stdio.h>

#include "dataset.h"
#include "util.h"

#define NUM_SIZE 7

void annotate_targets(dataset_t *tds, char* annotation_file);

#endif
