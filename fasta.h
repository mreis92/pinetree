#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>

#include "dataset.h"
#include "util.h"
#include "constants.h"

dataset_t *parse_fasta(char *filename);

#endif
