#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <stdio.h>

#include "constants.h"
#include "dataset.h"
#include "strmap.h"
#include "util.h"

void annotate_targets(dataset_t *tds, char* annotation_file);

#endif
