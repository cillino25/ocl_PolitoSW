#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <CL/cl.h>

#include "utils.h"



int getKernelMax(BLOCK * kernels, int size){
	int i, tmp = 0;
	for (i = 0; i < size; i++){
		//printf("Kernel %d: %d\n", i, kernels[i].h_max);			
		if (tmp < kernels[i].h_max)
			tmp = kernels[i].h_max;
	}
	return tmp;
}

int getKernelMaxLinear(BLOCK_LINEAR_SHORT * kernels, int size){
	int i, tmp = 0;

	for (i = 0; i < size; i++){
		//printf("Kernel %d: %d\n", i, kernels[i].h_max);			
		if (tmp < kernels[i].h_max)
			tmp = kernels[i].h_max;
	}
	return tmp;
}

// Does not work :-(
int getKernelMaxType(void * kernels, ALGORITHM type, int size){
	int i, tmp = 0;

	switch (type){
	case LINEAR:
		for (i = 0; i < size; i++){
			//printf("Kernel %d: %d\n", i, kernels[i].h_max);			
			if (tmp < ((BLOCK_LINEAR_SHORT*)kernels)[i].h_max)
				tmp = ((BLOCK_LINEAR_SHORT*)kernels)[i].h_max;
			break;
	default:
		break;
		}
	}
	return tmp;
}


int *get_buffer(char *sequence, char *symbols, int worker, int *buffer_size){
	int *buffer;
	char element;
	int value, i, j;
	int symbols_number = strlen(symbols);
	int size = strlen(sequence) + 1;
	if ((size%worker) != 0)
		size += worker - size%worker;	
	buffer = (int*)malloc(sizeof(int)*size);
	buffer[0] = 0;
	for (i = 1; i < size; i++){
		element = sequence[i - 1];
		value = symbols_number;
		for (j = 0; j < symbols_number && value == symbols_number; j++){
			if (symbols[j] == element)
				value = j;
		}
		buffer[i] = value+1;
	}
	*buffer_size = size;
	return buffer;
}

char *getCharBuffer(char *sequence, int sequence_length, char *symbols, int worker, size_t *buffer_size){
	char * buffer = NULL;
	char element = 0;
	int value = 0, i = 0, j = 0;
	int symbols_number = strlen(symbols);
	//size_t size = strlen(sequence) + 1;
	size_t size = (size_t) sequence_length + 1;

	/*if(sequence == NULL){
		*buffer_size = 0;
		return NULL;
	}*/

	if ((size%worker) != 0)
		size += worker - size%worker;
	
	buffer = (char*)malloc(sizeof(char)*size + 1);
	if(buffer == NULL){
		printf("Error allocating memory\nExiting..\n");
		return NULL;
	}
	buffer[0] = 5;
	for (i = 1; (i < size)&&(i < sequence_length+1); i++){
		element = sequence[i - 1];
		value = symbols_number;
		for (j = 0; j < symbols_number && value == symbols_number; j++){
			if (symbols[j] == element)
				value = j;
		}
		buffer[i] = (char)(value + 1);
	}
	
	for(; i < size; i++){
		buffer[i] = 5;
	}
	
	buffer[i] = 5;  // working
	//buffer[i] = '\0';
	
	*buffer_size = size;
	return buffer;
}

char * reallocGlobalBuffer(char *first, size_t first_size, char *append, size_t append_size, size_t *final_dest_size, unsigned long *offset){
	int i=0;
	
	//printf("DEBUG: reallocGlobalBuffer called with:\ninitial_dest_size=%ld\ndest  =", first_size);
	//for(i=0; i < first_size; i++) printf("%d", first[i]);
	//printf("\nsource_size=%ld\nsource=", append_size);
	//for(i=0; i < append_size; i++) printf("%d", append[i]);
	//printf("\n\n");
	
	//getchar();
	
	int final = first_size + append_size;
	char * tmp = (char*) calloc(final + 1, sizeof(char));
	//for(i=0; i<final+1; i++) tmp[i] = 5;
	*offset = (unsigned long)first_size;

	if(first == NULL){
		printf("reallocGlobalBuffer: dest pointer is NULL\n");
	}
	if(tmp == NULL){
		printf("Error allocating memory.\nExiting..\n");
		return NULL;
	}
	//printf("-------\ninitial_dest_size=%d, final_size=%d, source_size=%d\n", initial_dest_size, initial_dest_size+source_size, source_size);
	for(i = 0; i < first_size; i++){
		//printf("i=%2d: tmp[%d]=%d\n", i, i, dest[i]);
		tmp[i] = first[i];
	}
	
	//free(dest);

	//strncpy(tmp, source, initial_dest_size);
	
	
	
	for(i = 0; i < append_size; i++){
		tmp[i + first_size] = append[i];
		//printf("--dbg: i=%d: tmp[%d]=%d\n", i, i+initial_dest_size, tmp[i+initial_dest_size]);
	}
	
	//tmp[i+first_size] = '\0';   // working
	
	free(first);
	*final_dest_size = final;
	return(tmp);
}




BLOCK *get_block(int size, int pen){
	int i;
	BLOCK * tmp = (BLOCK*) malloc(sizeof(BLOCK) * size);

	for (i = 0; i < size; i++){
		tmp[i].up = 0;
		tmp[i].left = 0;
		tmp[i].diagonal = 0;
		tmp[i].h_max = 0;
		tmp[i].e = 0;
		tmp[i].f = 0;
		tmp[i].e_diagonal = 0;
		tmp[i].f_diagonal = 0;
		tmp[i].e_aux = 0;
		tmp[i].f_aux = 0;
		tmp[i].gQry = pen;
		tmp[i].gSbj = pen;
	}
	return tmp;
}

BLOCK_LINEAR_SHORT *get_block_linear(int size){
	int i;
	BLOCK_LINEAR_SHORT * tmp = (BLOCK_LINEAR_SHORT*)malloc(sizeof(BLOCK_LINEAR_SHORT)*size);

	for(i=0; i<size; i++){
		tmp[i].up = 0;
		tmp[i].left = 0;
		tmp[i].diagonal = 0;
		tmp[i].h_max = 0;
	}
	return tmp;
}


BLOCK_AUX *get_buffer_aux(int size){
	int i;
	BLOCK_AUX *buffer;
	buffer = (BLOCK_AUX*)malloc(sizeof(BLOCK_AUX)*size);
	for (i = 0; i < size; i++) {
		buffer[i].h = 0;
		buffer[i].e = 0;
		buffer[i].f = 0;
		buffer[i].gQry = 0;
		buffer[i].gSbj = 0;
		buffer[i].e_diagonal = 0;
		buffer[i].f_diagonal = 0;
		buffer[i].e_aux = 0;
		buffer[i].f_aux = 0;
	}
	return buffer;
}

BLOCK_AUX_AFFINE *get_buffer_aux_affine(int size){
	int i;
	BLOCK_AUX_AFFINE *buffer;
	buffer = (BLOCK_AUX_AFFINE*)malloc(sizeof(BLOCK_AUX_AFFINE)*size);
	for (i = 0; i < size; i++) {
		buffer[i].h = 0;
		buffer[i].e = 0;
		buffer[i].f = 0;
		buffer[i].e_diagonal = 0;
		buffer[i].f_diagonal = 0;
		buffer[i].e_aux = 0;
		buffer[i].f_aux = 0;
	}
	return buffer;
}

BLOCK_AFFINE *get_buffer_affine(int size){
	int i;
	BLOCK_AFFINE * tmp = (BLOCK_AFFINE*)malloc(sizeof(BLOCK_AFFINE)*size);

	for (i = 0; i<size; i++){
		tmp[i].up = 0;
		tmp[i].left = 0;
		tmp[i].diagonal = 0;
		tmp[i].h_max = 0;
		tmp[i].e = 0;
		tmp[i].f = 0;
		tmp[i].e_diagonal = 0;
		tmp[i].f_diagonal = 0;
		tmp[i].e_aux = 0;
		tmp[i].f_aux = 0;
	}
	return tmp;
}

BLOCK_AUX_LINEAR_SHORT *get_buffer_aux_linear(int size){
	int i;
	BLOCK_AUX_LINEAR_SHORT *buffer;
	buffer = (BLOCK_AUX_LINEAR_SHORT*)malloc(sizeof(BLOCK_AUX_LINEAR_SHORT)*size);
	for (i = 0; i < size; i++) {
		buffer[i].h = 0;
	}
	return buffer;
}

// Working
void divideEvenly2(int string_lenght, short segments, short *segment_length, short *complete_copiers, short *pad){
	short n, m, copiers;
	double tmp;

	if ((string_lenght % segments) == 0){
		n = string_lenght / segments;
		m = 0;
		copiers = segments;
	}
	else{
		copiers = MIN(segments - 1, string_lenght);
		tmp = (double)string_lenght / copiers;

		// Number of query elements copied by the single WI
		n = CEILING(tmp);

		// Needed for padding: last element will copy fewer elements than other WI
		m = string_lenght - n * copiers;

		if (m < 0)	// <==> WI are more than query elements
			m = 0;
	}

	*segment_length = n;
	*complete_copiers = copiers;
	*pad = m;
}

// 'complete_copiers' WI will copy 'elements' items, then the last element will copy 'pad_elements' item each
void divideEvenly(int string_lenght, short workers, short *complete_copiers, short *elements, short *pad){
	short copiers;
	short elems;
	short padding;

	if((string_lenght % workers) == 0){
		copiers = workers;
		elems = string_lenght / workers;
		padding = 0;
	}else if(string_lenght < workers){
		copiers = 0;
		elems = 0;
		padding = string_lenght;
	}else{
		copiers = MIN(workers - 1, string_lenght);			// copiers will copy elems items
		elems = CEILING((double) string_lenght / copiers) - 1;	// items each copier will copy
		padding = string_lenght - copiers * elems;
	}

	*complete_copiers = copiers;
	*elements = elems;
	*pad = padding;

}


// string_length must be divisible by workers
// Divides the string with best effort: with the maximum number of 16-wide blocks, plus a block padding, plus a single padding.

void divideEvenlyBlocks(int string_lenght, short workers, short blocks_width, short *complete_blocks, short *pad_block_copiers, short *single_pad_copiers){
	int remaining = 0;
	int loc_complete_blocks = 0;
	int loc_pad_block_copiers = 0;
	int loc_single_pad = 0;

	loc_complete_blocks = string_lenght / (workers * blocks_width);
	remaining = string_lenght - loc_complete_blocks * workers * blocks_width;
	
	
	
	loc_pad_block_copiers = remaining / 16;
	remaining = remaining - loc_pad_block_copiers * 16;
	
	// the rest will be divided among the other workers
	loc_single_pad = remaining / (workers - loc_pad_block_copiers);
	
	
	*complete_blocks = (short)loc_complete_blocks;		// how many times ALL the workers can copy a 16-wide value
	*pad_block_copiers = (short)loc_pad_block_copiers;		// how many workers can copy then a 16-wide value (surely < workers)
	*single_pad_copiers = (short)loc_single_pad;
}

void charv2intv(char *vchar, int *vint, int size){
	int i;
	for (i = 0; i < size; i++)
		vint[i] = (int)vchar[i];
	return;
}

void intv2charv(char *vchar, int *vint, int size){
	int i;
	for (i = 0; i < size; i++)
		vchar[i] = (char)vint[i];
	return;
}
