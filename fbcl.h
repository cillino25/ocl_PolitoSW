#ifndef _FBCL_H
#define _FBCL_H

#include <CL/cl.h>

#define MEM_SIZE 128
#define MAX_SOURCE_SIZE 1024

#define CL_MEM_TYPE_LOCAL	1<<3

typedef struct _fbcl_param FBCL_PARAM;

int fbclGetFirstPlatform(cl_platform_id *platform);
int fbclGetFirstDevice(cl_device_id *device, cl_device_type type);
int fbclGetFirstGPU(cl_device_id *device);
int fbclGetFirstCPU(cl_device_id *device);

// Gets device local memory in bytes
void fbclGetDeviceLocalMemory(cl_device_id device, cl_ulong * local_memory_size, cl_int * error);

void fbclGetMaxWorkGroupSize(cl_device_id device, size_t *size, cl_int *error, char * message);


char *fbclGetSourceFromFile(const char *file_name, size_t *source_size);



void fbclInitKernel(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, char * source, size_t source_size, cl_program *_program, cl_kernel *_kernel, char * function, int verbose);

/* 
:params: stringa che contiene i caratteri w [write] r [read] a [all] e descrive i parametri passati che sono nella forma [size_t dim, void* param] 
i parametri r -> sono zone di memoria in sola lettura
i parametri w -> sono zone di memoria in sola scrittura
i parametri a -> sono zone di memoria in lettura e scrittura
*/

// Original function + verbose
void fbclRun(cl_device_id device, int work_items, char *source, size_t source_size, char *function, int verbose, cl_double *execution_time, char *params_type, ...);

// Executes multiple alignments sequentially
// Provides support for LOCAL arguments as well
// TEMP version
int fbclRunSequential(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_program *_program, cl_kernel *_kernel, int work_items, char *source, size_t source_size, char *function, int verbose, cl_double *execution_time, char *params_type, ...);

// Executes multiple alignments in parallel. In particular, this function aligns many queries with one subject (TO DO: many queries with many sujects, efficiently managed)
// Provides support for LOCAL arguments as well
int fbclRunParallelSingle(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_kernel *_kernel, int work_items, char *function, int verbose, cl_double *execution_time, char *params_type, ...);

int fbclRunParallel(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_kernel *_kernel, int work_groups, int work_items, char *function, int verbose, cl_double *execution_time, char *params_type, ...);


// Very very original
//void fbclRun(cl_device_id device, int work_items, char *source, size_t source_size, char *function, double *execution_time, char *params_type, ...);

// TEST
// Executes runs instances of the kernel on the SAME command queue
// NB: executes ON THE SAME BUFFERS as well, so each kernel will override the previous kernel results
// (since the same output buffer is used)
void fbclSequentialRunSameQueue(cl_device_id device, int work_items, char *source, size_t source_size, char *function, int runs, int verbose, cl_double *execution_time, char *params_type, ...);



//fbclDestroy(&_context, &_program, &_kernel, &_queue);
int fbclDestroy(cl_context *_context, cl_program *_program, cl_kernel *_kernel, cl_command_queue *_queue);

cl_uint fbclNumberOfPlatforms();
cl_uint fbclNumberOfDevices(cl_platform_id platform, cl_device_type type);

void fbclDebug();

void fbclPrintDeviceInfo(cl_device_id device);

void fbclErrorTranslate(cl_int value, char *message);
void fbclProgramBuildInfo(cl_program program, cl_device_id device);

// Returns statistics on the given parameters
void fbclVariadicGeneral(long int *somma, int *massimo, int *minimo, int count, ...);

// Same as above, but pads the return values to 16 (in order to minimize memory accesses)
void fbclVariadicGeneral16(long int *somma, int *massimo, int *minimo, int count, ...);

// Same as above, but pads the return values to 32 (in order to fullfill an entire compute unit)
void fbclVariadicGeneral32(long int *somma, int *massimo, int *minimo, int count, ...);



#endif
