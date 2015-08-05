#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "file_handler.h"


//const char * query_file = "queries.txt";
//const char * subject_file = "subjects.txt";

/*
void main(){
	char ** sequences;
	int size, i, max;
	printf("Reading file %s\n", query_file);
	
	int l=0;
	sequences = readSequencesFromFile(query_file, &l, &size, &max);
	
	printf("Total: lines=%d, size=%d, max length=%d\n", l, size, max);
	for(i=0; i<l; i++){
		printf("line %d: %s\n", i, sequences[i]);
	}
	
	printf("\n\nIndexed version: Subjects file\n");
	SEQUENCE * s = readIndexedSequencesFromFile(subject_file , &l, &size, &max);
	printf("Total: lines=%d, size=%d, max length=%d\n", l, size, max);
	for(i=0; i<l; i++){
		printf("line %d: index=%d, length=%d, sequence=%s\n", i, s[i].index, s[i].length, s[i].seq);
	}
	
	printf("Sorting sequence by length:\n");
	sortSequencesByLength(s, l, 'd');
	for(i=0; i<l; i++){
		printf("line %d: index=%d, length=%d, sequence=%s\n", i, s[i].index, s[i].length, s[i].seq);
	}
	
}*/

char ** readSequencesFromFile(const char * filename, int * lines, int * file_size, int * maxLength){
	int i,j;
	int chars=0;
	int l=0;
	int size=0;
	int max=0;
	char c=0;
	FILE * fp;
	char buf[MAX_LINE];
	
	
	char ** seq = (char**)malloc(MAX_LINE * sizeof(char *));
	for(i=0; i<MAX_SEQS; i++){
		seq[i] = (char*)malloc(MAX_LINE * sizeof(char));
	}
	
	if((fp = fopen(filename, "r")) == NULL){
		printf("Error opening file %s\n", filename);
		return NULL;
	}
	
	while((( c = getc(fp)) != EOF )&&( chars < MAX_LINE )){
		
		if(c != '\n'){
			
			buf[chars] = c;
			chars++;
			
		}else{
			
			if(chars==0)
				continue;
			
			strncpy(seq[l], buf, chars);
			seq[l][chars] = '\0';
			
			if(max < chars) max = chars;
			size += chars;
			l++;
			chars = 0;
		}
	}
	fclose(fp);
	
	*lines = l;
	*file_size = size;
	*maxLength = max;
	
	char ** tmp = (char**) malloc(max * sizeof(char*));
	for(i=0; i<l; i++){
		tmp[i] = (char*) malloc(max * sizeof(char));
		
		for(j=0; j<max; j++){
			tmp[i][j] = seq[i][j];
		}
	}
	
	
	
	return(tmp);
}


SEQUENCE * readIndexedSequencesFromFile(const char * filename, int * lines, int * file_size, int * maxLength){
	int i=0,j=0;
	int chars=0;
	int l=0;
	int size=0;
	int max=0;
	char c=0;
	FILE * fp = NULL;
	char * buf = (char*) calloc(MAX_LINE, sizeof(char));
	
	
	SEQUENCE * seq = (SEQUENCE*)calloc(MAX_LINE, sizeof(SEQUENCE));
	for(i=0; i<MAX_SEQS; i++){
		seq[i].seq = (char*)calloc(MAX_LINE, sizeof(char));
		seq[i].index = 0;
		seq[i].seq_length = 0;
	}
	
	if((fp = fopen(filename, "r")) == NULL){
		printf("Error opening file %s\n", filename);
		return NULL;
	}
	
	while((( c = getc(fp)) != EOF )&&( chars < MAX_LINE )){
		
		
		if(c != '\n' && c != '\r'){
			if(chars == 0){
				char idx[MAX_INDEX_CHARS];
				int k=0;
				
				do{
					idx[k++] = c;
					//printf("DEBUG: idx[%d] = %c\n", k-1, idx[k-1]);
				}while((((c = getc(fp)) != ' ')&&(c != '\t'))&&( k<MAX_INDEX_CHARS ));
				
				if(k == MAX_INDEX_CHARS){
					printf("** WARNING: on %s : initial index of line %d exceeding %d chars (i.e. %d): maybe something's wrong.. **\n", filename, l, MAX_INDEX_CHARS, k);
					printf("\n");
					for(i=0; i<k; i++) printf("%c", idx[i]);
				}
				
				seq[l].index = atoi(idx);
				
				for(i=0; i<MAX_INDEX_CHARS; i++) idx[i] = '\0';
				//printf("DEBUG: line %d --> index %d\n", l, seq[l].index);
				c = getc(fp);
				
			}
			
			buf[chars] = c;
			chars++;
			
		}else{
			
			if(chars==0)
				continue;
			
			seq[l].seq_length = chars;
			strncpy(seq[l].seq, buf, chars);
			seq[l].seq[chars] = '\0';
			
			
			
			if(max < (chars-1)) max = chars;
			size += chars;
			l++;
			chars = 0;
		}
	}
	fclose(fp);
	
	if(chars == MAX_LINE){
		printf("** WARNING: line %d exceeded maximum dimensions. ***\n", l);
	}
	
	*lines = l;
	*file_size = size;
	*maxLength = max;
	
	SEQUENCE * tmp = (SEQUENCE *) calloc(l, sizeof(SEQUENCE));

	for(i=0; i<l; i++){
		tmp[i].seq = (char*) malloc(max * sizeof(char) + 1);
		strncpy(tmp[i].seq, seq[i].seq, max);
		tmp[i].seq[max] = '\0';
		
		tmp[i].index = seq[i].index;
		tmp[i].seq_length = seq[i].seq_length;
	}
	
	free(buf);
	for(i=0; i<MAX_SEQS; i++) free(seq[i].seq);
	free(seq);
	
	return(tmp);
}


void printSequence(SEQUENCE * s, int entries){
	int i=0, j=0;
	for(i=0; i < entries; i++){
		printf("\n* Index: %d\n* Length: %d\n* Sequence: ", s[i].index, s[i].seq_length);
		for(j=0; j<s[i].seq_length; j++)
			printf("%c", s[i].seq[j]);
		printf("\n\n");
	}
}


int seqCompareDescending(const void * a, const void * b){
	return ((*(SEQUENCE *)b).seq_length - (*(SEQUENCE *)a).seq_length);
}

int seqCompareAscending(const void * a, const void * b){
	return ((*(SEQUENCE *)a).seq_length - (*(SEQUENCE *)b).seq_length);
}

void sortSequencesByLength(SEQUENCE *s, int n, const char c){
	
	if(c == 'd')
		qsort(s, n, sizeof(SEQUENCE), seqCompareDescending);
	else if (c == 'a')
		qsort(s, n, sizeof(SEQUENCE), seqCompareAscending);
	else
		printf("** WARNING: Invalid sorting order requested. Not sorting.. **\n");
	
	return;
}