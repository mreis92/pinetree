#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "constants.h"
#include "types.h"

#define STRMATCH(s, n) strcmp(s, n) == 0 

void *safe_malloc(size_t);
void *safe_calloc(size_t, size_t);
void *safe_realloc(void *, size_t);
void safe_free(void *);
FILE *safe_fopen(const char *filename, const char *mode);
void safe_fclose(FILE *file);
FILE *safe_popen(const char *command, const char *mode);
void safe_pclose(FILE *file);
void safe_remove(char *filename);
char *get_system_time();
uint bin_code(char c);
uint int_code(char c);
char char_code(uint c);
char *reverse_complement(char *seq);
char *complement(char *seq);
char *reverse_sequence(char *seq);
char *get_string(FILE *, uint);
char *str_concat(char *, char *);
char *prepend_char(char *, char);
float fmax3(float, float, float);
int count_occurrences(char *string, char c);
void error(char *message);

#endif
