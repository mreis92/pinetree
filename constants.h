#ifndef CONSTANTS_H
#define CONSTANTS_H

#define VERSION 1.0

/*Program execution constants */
#define NUM_PROCESSORS 1
#define GAP -1
#define MATCH 1
#define MISMATCH -1
#define GU -0.5
#define SEED_START 2 
#define SEED_STOP 12
#define CR_START 9
#define CR_STOP 11
#define NUM_ERRORS 0
#define MARKOV_ORDER 3
#define EVALUE 40
#define C_THRESHOLD -1.1
#define A_THRESHOLD 2.3
#define E_THRESHOLD 7

#define PSEUDOCOUNTS 0.001

/* Default size for strings */
#define LONGBUF 1024
#define BUFSIZE 200
#define FASTASIZE 81
#define IDSIZE 10

/* Nucleotide parameters */
#define _A_ 0
#define _C_ 1
#define _T_ 2
#define _G_ 3
#define _N_ 4
#define NUCLEOTIDES 4

/* Constants for pretty print */
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

#endif
