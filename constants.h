#ifndef CONSTANTS_H
#define CONSTANTS_H

/*Program execution constants */
#define NUM_PROCESSORS 1
#define GAP -1
#define MATCH 1
#define MISMATCH -1
#define GU -0.5
#define SEED_PENALTY -0.5
#define SEED_START 2 
#define SEED_STOP 12
#define CR_START 9
#define CR_STOP 11
#define NUM_ERRORS 1
#define MARKOV_ORDER 6 
#define C_THRESHOLD 3.5
#define A_THRESHOLD 16
#define E_THRESHOLD 7

/* If you want to allow N errors in your sequences, define NUM_ERRORS as N+1 */
#define NUM_ERRORS 1

/* Default size for strings */
#define BUFSIZE 200
#define FASTASIZE 81

/* Nucleotide parameters */
#define _A_ 0
#define _C_ 1
#define _T_ 2
#define _G_ 3
#define _N_ 4
#define NUCLEOTIDES 4

#endif