#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "RNAplfold.h"

int unpaired;

/*--------------------------------------------------------------------------*/
double RNAplfold(int argc, char *argv[], char* id, char* seq, int region, int flanks){
	struct        RNAplfold_args_info args_info;
	char          fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH], *c, *structure, *ParamFile, *ns_bases, *rec_sequence, *rec_id, *orig_sequence;
	int           i, length, sym, winsize, pairdist;
	float         cutoff;
	int           tempwin, temppair, tempunpaired;
	FILE          *pUfp = NULL, *spup = NULL;
	double        **pup = NULL; /*prob of being unpaired, lengthwise*/
	int           noconv, plexoutput, simply_putout, openenergies, binaries;
	plist         *pl, *dpp = NULL;
	double        betaScale;
	double        score;
	pf_paramT     *pf_parameters;
	model_detailsT  md;

	dangles       = 2;
	cutoff        = 0.01;
	winsize       = 70;
	pairdist      = 0;
	unpaired      = 0;
	betaScale     = 1.;
	simply_putout = plexoutput = openenergies = noconv = 0;binaries=0;
	tempwin       = temppair = tempunpaired = 0;
	structure     = ParamFile = ns_bases = NULL;
	rec_id        = rec_sequence = orig_sequence = NULL;
	pf_parameters = NULL;

	set_model_details(&md);

	/*
	#############################################
	# check the command line parameters
	#############################################
	*/
	if(RNAplfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
	/* temperature */
	if(args_info.temp_given)              temperature = args_info.temp_arg;
	/* do not take special tetra loop energies into account */
	if(args_info.noTetra_given)           md.special_hp = tetra_loop=0;
	/* set dangle model */
	if(args_info.dangles_given){
		if((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
			warn_user("required dangle model not implemented, falling back to default dangles=2");
		else
			md.dangles = dangles = args_info.dangles_arg;
	}
	/* do not allow weak pairs */
	if(args_info.noLP_given)              md.noLP = noLonelyPairs = 1;
	/* do not allow wobble pairs (GU) */
	if(args_info.noGU_given)              md.noGU = noGU = 1;
	/* do not allow weak closing pairs (AU,GU) */
	if(args_info.noClosingGU_given)       md.noGUclosure = no_closingGU = 1;
	/* do not convert DNA nucleotide "T" to appropriate RNA "U" */
	if(args_info.noconv_given)            noconv = 1;
	/* set energy model */
	if(args_info.energyModel_given)       energy_set = args_info.energyModel_arg;
	/* take another energy parameter set */
	if(args_info.paramFile_given)         ParamFile = strdup(args_info.paramFile_arg);
	/* Allow other pairs in addition to the usual AU,GC,and GU pairs */
	if(args_info.nsp_given)               ns_bases = strdup(args_info.nsp_arg);
	/* set the maximum base pair span */
	if(args_info.span_given)              pairdist = args_info.span_arg;
	/* set the pair probability cutoff */
	if(args_info.cutoff_given)            cutoff = args_info.cutoff_arg;
	/* set the windowsize */
	if(args_info.winsize_given)           winsize = args_info.winsize_arg;
	/* set the length of unstructured region */
	if(args_info.ulength_given)           unpaired = args_info.ulength_arg;
	/* compute opening energies */
	if(args_info.opening_energies_given)  openenergies = 1;
	/* print output on the fly */
	if(args_info.print_onthefly_given)    simply_putout = 1;
	/* turn on RNAplex output */
	if(args_info.plex_output_given)       plexoutput = 1;
	/* turn on binary output*/
	if(args_info.binaries_given)          binaries = 1;

	if(args_info.betaScale_given)         betaScale = args_info.betaScale_arg;

	/* check for errorneous parameter options */
	if((pairdist < 0) || (cutoff < 0.) || (unpaired < 0) || (winsize < 0)){
		RNAplfold_cmdline_parser_print_help();
		exit(EXIT_FAILURE);
	}

	/* free allocated memory of command line data structure */
	RNAplfold_cmdline_parser_free(&args_info);

	/*
	#############################################
	# begin initializing
	#############################################
	*/
	if (ParamFile != NULL)
	read_parameter_file(ParamFile);

	if (ns_bases != NULL) {
		nonstandards = space(33);
		c=ns_bases;
		i=sym=0;
		if (*c=='-') {
			sym=1; c++;
		}
		while (*c!='\0') {
			if (*c!=',') {
				nonstandards[i++]=*c++;
				nonstandards[i++]=*c;
				if ((sym)&&(*c!=*(c-1))) {
					nonstandards[i++]=*c;
					nonstandards[i++]=*(c-1);
				}
			}
			c++;
		}
	}

	/* check parameter options again and reset to reasonable values if needed */
	if(openenergies && !unpaired) unpaired  = 31;
	if(pairdist == 0)             pairdist  = winsize;
	if(pairdist > winsize){
		fprintf(stderr, "pairdist (-L %d) should be <= winsize (-W %d);"
		"Setting pairdist=winsize\n",pairdist, winsize);
		pairdist = winsize;
	}
	if(dangles % 2){
		warn_user("using default dangles = 2");
		dangles = 2;
	}

	rec_id = id;
	rec_sequence = seq;


	/*
	########################################################
	# init everything according to the data we've read
	########################################################
	*/
	if(rec_id){
		(void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
	}
	else fname[0] = '\0';

	length    = (int)strlen(rec_sequence);
	structure = (char *) space((unsigned) length+1);

	/* convert DNA alphabet to RNA if not explicitely switched off */
	if(!noconv) str_DNA2RNA(rec_sequence);
	/* store case-unmodified sequence */
	orig_sequence = strdup(rec_sequence);
	/* convert sequence to uppercase letters only */
	str_uppercase(rec_sequence);

	/*
	########################################################
	# done with 'stdin' handling
	########################################################
	*/

	if(length > 1000000){
		if(!simply_putout && !unpaired){
			printf("Switched to simple output mode!!!\n");
			simply_putout = 1;
		}
	}
	if(unpaired && simply_putout){
		printf("Output simplification not possible if unpaired is switched on\n");
		simply_putout = 0;
	}

	/* restore winsize if altered before */
	if(tempwin != 0){
		winsize = tempwin;
		tempwin = 0;
	}
	/* restore pairdist if altered before */
	if(temppair != 0){
		pairdist = temppair;
		temppair = 0;
	}
	/* restore ulength if altered before */
	if(tempunpaired != 0){
		unpaired      = tempunpaired;
		tempunpaired  = 0;
	}

	/* adjust winsize, pairdist and ulength if necessary */
	if(length < winsize){
		/*fprintf(stderr, "WARN: window size %d larger than sequence length %d\n", winsize, length);*/
		tempwin = winsize;
		winsize = length;
		if (pairdist>winsize) {
			temppair=pairdist;
			pairdist=winsize;
		}
		if (unpaired>winsize) {
			tempunpaired=unpaired;
			unpaired=winsize;
		}
	}

	/*
	########################################################
	# begin actual computations
	########################################################
	*/

	if (length >= 5){
		/* construct output file names */
		char fname1[FILENAME_MAX_LENGTH], fname2[FILENAME_MAX_LENGTH], fname3[FILENAME_MAX_LENGTH], fname4[FILENAME_MAX_LENGTH], fname_t[FILENAME_MAX_LENGTH];

		strcpy(fname_t, (fname[0] != '\0') ? fname : "plfold");

		strcpy(fname1, fname_t);
		strcpy(fname2, fname_t);
		strcpy(fname3, fname_t);
		strcpy(fname4, fname_t);
		strcpy(ffname, fname_t);

		strcat(fname1, "_lunp");
		strcat(fname2, "_basepairs");
		strcat(fname3, "_uplex");
		if(binaries){
			strcat(fname4, "_openen_bin");
		}
		else{
			strcat(fname4, "_openen");
		}
		strcat(ffname, "_dp.ps");

		pf_parameters = get_boltzmann_factors(temperature, betaScale, md, -1);

		if(unpaired > 0){
			pup       =(double **)  space((length+1)*sizeof(double *));
			pup[0]    =(double *)   space(sizeof(double)); /*I only need entry 0*/
			pup[0][0] = unpaired;
		}

		pUfp = spup = NULL;
		if(simply_putout){
			spup = fopen(fname2, "w");
			pUfp = (unpaired > 0) ? fopen(fname1, "w") : NULL;

			pl = pfl_fold_par(rec_sequence, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup, pf_parameters);

			if(pUfp != NULL)  fclose(pUfp);
			if(spup != NULL)  fclose(spup);
		}
		else{
			pl = pfl_fold_par(rec_sequence, winsize, pairdist, cutoff, pup, &dpp, pUfp, spup, pf_parameters);
			/*PS_dot_plot_turn(orig_sequence, pl, ffname, pairdist);*/
			if (unpaired > 0){
				if(plexoutput){
					pUfp = fopen(fname3, "w");
					putoutphakim_u(pup,length, pUfp);
					fclose(pUfp);
				}
				/*pUfp = fopen(openenergies ? fname4 : fname1, "w");*/
				if(binaries){
					putoutpU_prob_bin_par(pup, length, unpaired, pUfp, openenergies, pf_parameters);
				}
				else{
					/*putoutpU_prob_par(pup, length, unpaired, pUfp, openenergies, pf_parameters);*/
				}
				/*fclose(pUfp);*/
			}
		}

		score = -log(pup[region+flanks][region])*pf_parameters->kT/1000;
		
		for (i=1; i<=length; i++){
	    		free(pup[i]);
	  	}

		if(pl) free(pl);
		if(unpaired > 0){
			free(pup[0]);
			free(pup);
		}

		free(pf_parameters);
	}
	(void) fflush(stdout);

	/* clean up */
	if(rec_id) free(rec_id);
	free(rec_sequence);
	free(orig_sequence);
	free(structure);
	rec_id = rec_sequence = orig_sequence = NULL;

	return score;
}

void putoutphakim_u(double **pU,int length, FILE *fp) {
	/*put out Fopen in dekacalories per mol, and F(cond,open) also in dekacal*/
	int k;

	float RT = (temperature+K0)*GASCONST;
	float p0;
	float pdep;
	int f0;
	int fdep;

	fprintf(fp,"#energy necessary to unpair as well as to unpair if i-1 is unpaired also, if i+1 is unpaired also in dekacal/mol\n");
	for (k=1; k<=length; k++){
		fprintf(fp,"%d\t",k);
		p0=pU[k][1];
		f0=(int) -RT*log(p0)/10;
		fprintf(fp,"%d\t", f0);
		if (k>1) {
			pdep=pU[k][2]/pU[k-1][1];
			fdep=(int) -RT*log(pdep)/10;
			fprintf(fp,"%d\t",fdep);
		}
		else  fprintf(fp,"-0\t");
			if (k<length) {
			pdep=pU[k+1][2]/pU[k+1][1];
			fdep=(int) -RT*log(pdep)/10;
			fprintf(fp,"%d\t",fdep);
		}
		else  fprintf(fp,"-0\t");
		fprintf(fp,"\n");
	}

		fflush(fp);
}

