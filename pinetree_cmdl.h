#include <stdio.h>

#include "align.h"
#include "pinetree_utils.h"
#include "util.h"

void initialize_parameters(char* filename, pinetree_args* args);
void common_usage();
void print_usage();
pinetree_args* read_cml_arguments(int argc, char **argv);