#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include "fasta.h"

/* Parses a FASTA file and returns a dataset from it.
	The headers will be the id for each sequence */
dataset_t *parse_fasta(char *filename)
{
	FILE *fasta = safe_fopen(filename, "r");
	
	int seqn = 0, seenheader = 0;
	char **sequences = NULL, **ids = NULL;
	char *buffer = NULL, *temp_buffer = NULL;

	while (!feof(fasta)) {
		char *current_str = get_string(fasta, FASTASIZE);
		switch (current_str[0]) {
		case '>':
			if (seenheader) {
				sequences[seqn++] = buffer;
				buffer = NULL;
			}

			ids = safe_realloc(ids, sizeof(char **) * (seqn + 1));
			sequences =
			    safe_realloc(sequences,
					 sizeof(char **) * (seqn + 1));

			/* removing the '>' character from the beginning of the string */
			ids[seqn] = strdup(current_str + 1);
			seenheader = 1;
			break;
		default:

			if (!seenheader)
				error("Invalid FASTA file.\n");

			temp_buffer = str_concat((buffer == NULL) ? "" : buffer, current_str);
			safe_free(buffer);
			buffer = temp_buffer;
			break;
		}

		safe_free(current_str);
	}

	sequences[seqn++] = buffer;
	fclose(fasta);
	
	return create_dataset(seqn, sequences, ids);
}
