/* Last changed Time-stamp: <2008-07-02 17:10:04 berni> */
/*
                  Ineractive Access to folding Routines

                  c Ivo L Hofacker
                  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "PS_dot.h"
#include "energy_const.h"
#include "read_epars.h"
#include "LPfold.h"
#include "params.h"
#include "RNAplfold_cmdl.h"

void putoutphakim_u(double **pU,int length, FILE *fp);

/*--------------------------------------------------------------------------*/
double RNAplfold(int argc, char *argv[], char *id, char *seq, int region, int flanks);
