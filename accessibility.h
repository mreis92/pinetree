#ifndef ACCESSIBILITY_H
#define ACCESSIBILITY_H

#include  <stdio.h>

#include "util.h"
#include "RNAup.h"
#include "RNAplfold.h"

#define UPSTREAM 17
#define DOWNSTREAM 13
#define FLANK 100
#define HEADER 1

double calculate_accessibility(char* mode, char *gene_seq, int start, int align_length);
double accessibility_RNAup(char *gene, int region, int end_pos);
double accessibility_RNAplfold(char *seq, int region, int flanks);
/*double accessibility_RNAhybrid(char* gene, char* mirna);*/

#endif

