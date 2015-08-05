
#ifndef _UTILS_H
#define	_UTILS_H

#define BUFFER_SIZE 50
#define KERNEL_RUNS 10

#define YES 1
#define NO 0

#define WINDOWS	1
#define WINDOWS_PROFILE	1

#ifdef WINDOWS_PROFILE
#define PROFILE_START()\
    QueryPerformanceCounter(&start);
#define PROFILE_STOP()\
    QueryPerformanceCounter(&end);\
	exec_time = (double)(end.QuadPart - start.QuadPart) * 1000 / frequency.QuadPart;

#define TOT_PROFILE_START()\
    QueryPerformanceCounter(&tot_start);
#define TOT_PROFILE_STOP()\
    QueryPerformanceCounter(&tot_end);\
	tot_exec_time = (double)(tot_end.QuadPart - tot_start.QuadPart) * 1000 / frequency.QuadPart;

#else

#define PROFILE_START()\
    gettimeofday(&tv_start, NULL);
#define PROFILE_STOP()\
    gettimeofday(&tv_end, NULL);\
    exec_time = (double)(tv_end.tv_usec - tv_start.tv_usec)/1000 + (double)(tv_end.tv_sec - tv_start.tv_sec)*1000;

#define TOT_PROFILE_START()\
    gettimeofday(&tv_tot_start, NULL);
#define TOT_PROFILE_STOP()\
    gettimeofday(&tv_tot_end, NULL);\
    tot_exec_time = (double)(tv_tot_end.tv_usec - tv_tot_start.tv_usec)/1000 + (double)(tv_tot_end.tv_sec - tv_tot_start.tv_sec)*1000;
#endif


// Empirically defined, measures the number of couples (q,s) that will be loaded into one work-group at a time,
// in order to minimize memory transfers (to fill device memory and allow it to do as much computation as possible).
// This number heavily depends on the size of query and subjects, so it could be determined by profiling the applications
// or statically evaluating the query and subject databases.
#define	RUNS_PER_KERNEL		10

// theorically BEST number, since in TITAN we have 32 threads/warp (warp: indipendend flow of control,
// each warp is independent on the others when control flow divides - eg. when an IF is encountered)
#define WORKERS_PER_GROUP	128

// TITAN specific informations:
#define CUDA_CORES								2688
#define THREAD_PER_WARP							32				// 32 threads all executing the same instructions, on different data
#define WARPS_PER_MULTIPROCESSOR				64
#define THREADS_PER_MULTIPROCESSOR				2048
#define MAX_THREAD_BLOCKS_PER_MULTIPROCESSOR	16				// MAX WORK GROUPS!!!
#define MULTIPROCESSORS							14				// OpenCL Compute Units
#define GLOBAL_MEMORY							6442450944
#define LOCAL_MEMORY							49151


// Local work group size = 32,
// maximum number of work groups = 15 * 64 = 960		==> theoretically 960 (Q,S) parallel alignments!!
// OR maximum = 15 * 16 = 240 ???
//  
// OpenCL dimensions: 1024 * 1024 * 64
//
//		In one dimension we can map (1024 / 32) = 32 work groups.
//		==> In order to map 960 work groups we need (960 / 32) = 30 dimensions.

// OpenCL renaming for GTX Titan
#define TITAN_COMPUTE_UNITS						15
#define TITAN_MAX_WORK_GROUPS					240		// = 15 * 16
#define TITAN_LINEAR_DIMENSION					1024
#define TITAN_WG_PER_LINEAR_DIMENSION			32		// = 1024 / 32
#define TITAN_DIMENSIONS						8		//	240 / 32 = 7.5 ==> 8




// Intel HD 4600:
#define HD4600_COMPUTE_UNITS					20		// not really compute units, since local memory is shared between the both slices (ie. the 20 CUs)
														// Logic barriers allows up to 16 work-groups to be executed independently per subslice ?

#define HD4600_MAX_WORK_GROUPS					32		// 16 work-groups * 2 subslices

#define HD4600_LINEAR_DIMENSION					512
#define HD4600_WG_PER_LINEAR_DIMENSION			16		// = 1024 / 32 (32=workers_per_work-group)
#define HD4600_DIMENSIONS						2		// 32 / 16
#define HD4600_VERTICAL_GROUPS                  2


#define HD4600		1
//#define GTX_TITAN	1

#ifdef HD4600
#define COMPUTE_UNITS				HD4600_COMPUTE_UNITS
#define MAX_WORK_GROUPS				HD4600_MAX_WORK_GROUPS
#define	WG_PER_LINEAR_DIMENSION		HD4600_WG_PER_LINEAR_DIMENSION
#define DIMENSIONS					HD4600_DIMENSIONS
#define LINEAR_DIMENSION            HD4600_LINEAR_DIMENSION
#define VERTICAL_GROUPS             HD4600_VERTICAL_GROUPS
#endif

#ifdef GTX_TITAN
#define COMPUTE_UNITS				TITAN_COMPUTE_UNITS
#define MAX_WORK_GROUPS				TITAN_MAX_WORK_GROUPS
#define	WG_PER_LINEAR_DIMENSION		TITAN_WG_PER_LINEAR_DIMENSION
#define DIMENSIONS					TITAN_DIMENSIONS
#define LINEAR_DIMENSION            TITAN_LINEAR_DIMENSION
#endif





#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define CEILING_POS(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))
#define CEILING_NEG(X) ((X-(int)(X)) < 0 ? (int)(X-1) : (int)(X))
#define CEILING(X) ( ((X) > 0) ? CEILING_POS(X) : CEILING_NEG(X) )

/*************************************************************************************/
typedef enum{ LINEAR, AFFINE, DGS } ALGORITHM;

// DNA sub matrix, match= 2, unmatch= -1
const int dna[36] = {
	0, 0, 0, 0, 0, 0,
	0, 2, -1, -1, -1, 0,
	0, -1, 2, -1, -1, 0,
	0, -1, -1, 2, -1, 0,
	0, -1, -1, -1, 2, 0,
	0, 0, 0, 0, 0, 0
};



const char char_sub_matrix[36] = {//[45] = {
	0,  0,  0,  0,  0,  0,
	0,  2, -1, -1, -1,  0,
	0, -1,  2, -1, -1,  0,
	0, -1, -1,  2, -1,  0,
	0, -1, -1, -1,  2,  0,
	0,  0,  0 , 0,  0,  0//, 
	//0,0,0,0,0,0,0,0,0
};

typedef struct _info {
	int	reference_lenght;
	int	query_lenght;
	int symbol_number; 
	int block_size;
	int block_x;
	int block_y;	
	int penality1;
	int penality2;
} INFO;

typedef struct _global_info {
	cl_int symbol_number;
	cl_int block_size;
	cl_int penality1;
	cl_int penality2;
    
    cl_int groups_per_line;
    cl_int total_alignments;
    cl_int parallel_jobs;
} GLOBAL_INFO;


typedef struct _block_aux {
	int h;
	int e;
	int f;
	int e_diagonal;
	int f_diagonal;
	int e_aux;
	int f_aux;
	int gSbj;
	int gQry;
} BLOCK_AUX;

typedef struct _block {
	int up;
	int left;
	int diagonal;
	int h_max;
	int e;
	int f;
	int e_diagonal;
	int f_diagonal;
	int e_aux;
	int f_aux;
	int gSbj;
	int gQry;
} BLOCK;


/*
 *	Data structures for LINEAR algorithm
 */

typedef struct _block_aux_linear_short {
	cl_short h;
} BLOCK_AUX_LINEAR_SHORT;

typedef struct _block_linear_short {
	cl_short up;
	cl_short left;
	cl_short diagonal;
	cl_short h_max;
} BLOCK_LINEAR_SHORT;


typedef struct _control_short {
    cl_short ref_complete_blocks;
    cl_short ref_pad_block_copiers;
    cl_short ref_single_pad_copiers;
	cl_short qry_complete_blocks;
	cl_short qry_pad_block_copiers;
	cl_short qry_single_pad_copiers;
	cl_short sub_complete_copiers;
	cl_short sub_complete_elements;
	cl_short sub_pad_elements;
} CONTROL_SHORT;


typedef struct _control_parallel {
	cl_int	reference_length;
    cl_int  reference_padded_length;
	cl_int	query_length;
    cl_int  query_padded_length;
	cl_int  block_x;
	cl_int  block_y;

	cl_ulong ref_id;
	cl_ulong ref_offset;
	cl_int   ref_complete_blocks;
	cl_int   ref_pad_block_copiers;
	cl_int   ref_single_pad_copiers;

	cl_ulong qry_id;
	cl_ulong qry_offset;
	cl_int   qry_complete_blocks;
	cl_int   qry_pad_block_copiers;
	cl_int   qry_single_pad_copiers;

	cl_short sub_complete_copiers;
	cl_short sub_complete_elements;
	cl_short sub_pad_elements;
	
	// output
    cl_int gid_x;
    cl_int gid_y;
    cl_int iter;
    cl_int idx;
    
	cl_int score;
} CONTROL_PARALLEL;


/*
 *	Data structures for AFFINE algotihm
 */

typedef struct _block_aux_affine {
	cl_short h;
	cl_short e;
	cl_short f;
	cl_short e_diagonal;
	cl_short f_diagonal;
	cl_short e_aux;
	cl_short f_aux;
	/*cl_int h;
	cl_int e;
	cl_int f;
	cl_int e_diagonal;
	cl_int f_diagonal;
	cl_int e_aux;
	cl_int f_aux;*/
} BLOCK_AUX_AFFINE;

typedef struct _block_affine {
	cl_int up;
	cl_int left;
	cl_int diagonal;
	cl_int h_max;
	cl_int e;
	cl_int f;
	cl_int e_diagonal;
	cl_int f_diagonal;
	cl_int e_aux;
	cl_int f_aux;
} BLOCK_AFFINE;


/*
*	Data structures for DGS algotihm
*/

typedef struct _block_aux_dgs {
	cl_short h;
	cl_short e;
	cl_short f;
	cl_short e_diagonal;
	cl_short f_diagonal;
	cl_short e_aux;
	cl_short f_aux;
	cl_short gSbj;
	cl_short gQry;
	/*cl_int h;
	cl_int e;
	cl_int f;
	cl_int e_diagonal;
	cl_int f_diagonal;
	cl_int e_aux;
	cl_int f_aux;
	cl_short gSbj;
	cl_short gQry;*/
} BLOCK_AUX_DGS;

typedef struct _block_dgs {
	cl_int up;
	cl_int left;
	cl_int diagonal;
	cl_int h_max;
	cl_int e;
	cl_int f;
	cl_int e_diagonal;
	cl_int f_diagonal;
	cl_int e_aux;
	cl_int f_aux;
	cl_int gSbj;
	cl_int gQry;
} BLOCK_DGS;







void testSW(cl_device_id _device, int worker, ALGORITHM alg);
void originalTestSW(cl_device_id device, int worker, ALGORITHM alg);


// Sequential alignment of all the possibile couples (Q,S) made with the two files
// Content of out_file will be overwritten
void clAlignSWSequentialChar(cl_device_id device, int worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix);


// Single parallel execution
int clAlignSWParallelSingle(cl_device_id device, int max_worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix);

// Parallel alignment of all the possibile couples (Q,S) made with the two files
// Content of out_file will be overwritten.

// HP: subject db is MUCH BIGGER than query db!!!!
// If seen as a matrix of (s,q) couples, subjects are in the rows while queries in the columns, and the resulting matrix will be MUCH TALLER than wider
// For this reason, it is efficient to set the single work-group to scan the matrix VERTICALLY, that is with the same query but multiple subjects.
int clAlignSWParallelChar(cl_device_id device, int parallel_jobs, int max_worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix);



/************************************************** Util functions ********************************************************/

int getKernelMax(BLOCK * kernels, int size);
int getKernelMaxLinear(BLOCK_LINEAR_SHORT * kernels, int size);
int getKernelMaxType(void * kernels, ALGORITHM type, int size);


int *get_buffer(char *sequence, char *symbols, int worker, int *buffer_size);
char *getCharBuffer(char *sequence, int sequence_length, char *symbols, int worker, size_t *buffer_size);

// adds source to the end of dest, and sets the offset accordingly
char * reallocGlobalBuffer(char *first, size_t first_size, char *append, size_t append_size, size_t *final_dest_size, unsigned long *offset);


BLOCK *get_block(int size, int pen);
BLOCK_LINEAR_SHORT *get_block_linear(int size);

BLOCK_AUX *get_buffer_aux(int size);
BLOCK_AUX_LINEAR_SHORT *get_buffer_aux_linear(int size);

BLOCK_AFFINE *get_buffer_affine(int size);
BLOCK_AUX_AFFINE *get_buffer_aux_affine(int size);

void divideEvenly(int string_lenght, short workers, short *complete_copiers, short *elements, short *pad);		// Final version!
void divideEvenly2(int string_lenght, short segments, short *segment_length, short *complete_copiers, short *pad); // can be removed



void divideEvenlyBlocks(int string_lenght, short workers, short blocks_width, short *complete_blocks, short *pad_block_copiers, short *single_pad_copiers);

void charv2intv(char *vchar, int *vint, int size);
void intv2charv(char *vchar, int *vint, int size);


#endif	/* _UTILS_H */

