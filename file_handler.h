#ifndef _FILE_HANDLER_H
#define _FILE_HANDLER_H

#include <CL/cl.h>

//const char * query_file = "queries.txt";
//const char * subject_file = "subjects.txt";


#define MAX_LINE            50000    // max query/subject length     (columns)
#define MAX_SEQS            1000    // max number of lines per file (rows)
#define MAX_INDEX_CHARS     8       // maximum chars for the initial index

typedef struct sequence {
	char * seq;
	cl_int index;
	cl_int seq_length;
} SEQUENCE;

typedef struct _cl_sequence {
    char seq[MAX_SEQS];
    cl_int index;
    cl_int seq_length;
} CL_SEQUENCE;

/*
 *  Read from query file and writes them back into the char ** queries (array of strings).
 *  Returns a positive number representing the number of read queries on success,
 *  a non-positive (<= 0) number otherwise.
 *  Query and subject files are supposed to begin with an initial non-negative integer number
 *  representing the index of the sequence inside the file.
 * 
 *  NB: files must terminate with '\n'
 */

char ** readSequencesFromFile(const char * filename, int * lines, int * file_size, int * maxLength);

SEQUENCE * readIndexedSequencesFromFile(const char * filename, int * lines, int * file_size, int * maxLength);

void printSequence(SEQUENCE * s, int entries);

// Sorts the structure of sequences by length, in ascending (c='a') or descending order (c='d').
// It uses qsort() given by stdlib
void sortSequencesByLength(SEQUENCE *s, int n, const char c);

#endif