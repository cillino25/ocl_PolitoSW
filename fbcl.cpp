#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#include "fbcl.h"
#include "utils.h"

struct _fbcl_param{
	void *param;
	unsigned long param_size;
	cl_mem memory;
	cl_mem_flags memory_type;
	cl_mem * multi_mem;
};

// Original + verbose
void fbclRun(cl_device_id device, int work_items, char *source, size_t source_size, char *function, int verbose, double *execution_time, char *params_type, ...){
	va_list variadic_list;
	cl_context _context;
	cl_command_queue _queue;
	cl_mem _memory;
	cl_kernel _kernel;
	cl_program _program;
	cl_event _profiling_event;
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[1], local_work_size[1];
	
	//Pre analysis
	param_element = strlen(params_type);	
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	if(verbose){
		printf("Run kernel function: %s\nDevice:\n", function);
		fbclPrintDeviceInfo(device);
	}
	
	_context = clCreateContext(NULL, 1, &device, NULL, NULL, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("\nDEBUG: context\t%s\n", message);
	}
	
	_queue = clCreateCommandQueue(_context, device, CL_QUEUE_PROFILING_ENABLE, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: queue\t%s\n", message);		
	}
	
	_program = clCreateProgramWithSource(_context, 1, (const char **)&source, (const size_t *)&source_size, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: program\t%s\n", message);
	}

	//Build
	_ret = clBuildProgram(_program, 1, &device, NULL, NULL, NULL);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: build\t%s\n", message);
		fbclProgramBuildInfo(_program, device);
	}
	
	_kernel = clCreateKernel(_program, function, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: kernel\t%s\n", message);
	}
		
	//Memory INIT
	for (i = 0; i < param_element; i++){
		params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{
		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		if(verbose){
			printf("Param %i: Type %c Size %ld\n", i, params_type[i], params[i].param_size);
		}
		
		params[i].memory = clCreateBuffer(_context, params[i].memory_type, params[i].param_size, NULL, &_ret);
		
		if(verbose){
			fbclErrorTranslate(_ret, message);
			printf("DEBUG: memory\t%s\n", message);
		}
	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(_queue, params[i].memory, CL_FALSE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message); }
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		_ret = clSetKernelArg(_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: init\t%s\n", message); }
	}


	//Execute
	//_ret = clEnqueueTask(_queue, _kernel, 0, NULL, NULL);
	global_work_size[0] = work_items;
	local_work_size[0] = work_items;
	_ret = clEnqueueNDRangeKernel(_queue, _kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event);
	if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message); }

	clWaitForEvents(1, &_profiling_event);
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);

	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message); }
		}
	}

	//Destroy all
	_ret = clFlush(_queue);
	_ret = clFinish(_queue);
	_ret = clReleaseEvent(_profiling_event);
	_ret = clReleaseKernel(_kernel);
	_ret = clReleaseProgram(_program);	
	for (i = 0; i < param_element; i++){
		_ret = clReleaseMemObject(params[i].memory);
	}
	_ret = clReleaseCommandQueue(_queue);
	_ret = clReleaseContext(_context);
	free(params);

	if(verbose) printf("DEBUG: finish\n");
	return;
}
// Very original
/*
void fbclRun(cl_device_id device, int work_items, char *source, size_t source_size, char *function, double *execution_time, char *params_type, ...){
	va_list variadic_list;
	cl_context _context;
	cl_command_queue _queue;
	cl_mem _memory;
	cl_kernel _kernel;
	cl_program _program;
	cl_event _profiling_event;
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[1], local_work_size[1];
	
	//Pre analysis
	param_element = strlen(params_type);	
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	printf("Run kernel function: %s\nDevice:\n", function); fbclPrintDeviceInfo(device);
	_context = clCreateContext(NULL, 1, &device, NULL, NULL, &_ret);
	fbclErrorTranslate(_ret, message); printf("\nDEBUG: context\t%s\n", message);
	_queue = clCreateCommandQueue(_context, device, CL_QUEUE_PROFILING_ENABLE, &_ret);
	fbclErrorTranslate(_ret, message); printf("DEBUG: queue\t%s\n", message);		
	_program = clCreateProgramWithSource(_context, 1, (const char **)&source, (const size_t *)&source_size, &_ret);
	fbclErrorTranslate(_ret, message); printf("DEBUG: program\t%s\n", message);

	//Build
	_ret = clBuildProgram(_program, 1, &device, NULL, NULL, NULL);
	fbclErrorTranslate(_ret, message); printf("DEBUG: build\t%s\n", message); fbclProgramBuildInfo(_program, device);
	_kernel = clCreateKernel(_program, function, &_ret);
	fbclErrorTranslate(_ret, message); printf("DEBUG: kernel\t%s\n", message);
		
	//Memory INIT
	for (i = 0; i < param_element; i++){
		params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{
		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		printf("Param %i: Type %c Size %i\n", i, params_type[i], params[i].param_size);
		params[i].memory = clCreateBuffer(_context, params[i].memory_type, params[i].param_size, NULL, &_ret);
		fbclErrorTranslate(_ret, message); printf("DEBUG: memory\t%s\n", message);
	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(_queue, params[i].memory, CL_FALSE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message);
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		_ret = clSetKernelArg(_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		fbclErrorTranslate(_ret, message); printf("DEBUG: init\t%s\n", message);
	}


	//Execute
	//_ret = clEnqueueTask(_queue, _kernel, 0, NULL, NULL);
	global_work_size[0] = work_items;
	local_work_size[0] = work_items;
	_ret = clEnqueueNDRangeKernel(_queue, _kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event);
	fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message);

	clWaitForEvents(1, &_profiling_event);
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);

	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message);
		}
	}

	//Destroy all
	_ret = clFlush(_queue);
	_ret = clFinish(_queue);
	_ret = clReleaseKernel(_kernel);
	_ret = clReleaseProgram(_program);	
	for (i = 0; i < param_element; i++){
		_ret = clReleaseMemObject(params[i].memory);
	}
	_ret = clReleaseCommandQueue(_queue);
	_ret = clReleaseContext(_context);
	free(params);

	printf("DEBUG: finish\n");
	return;
}*/

void fbclSequentialRunSameQueue(
	cl_device_id device, int work_items, char *source, size_t source_size, char *function, int runs, int verbose, cl_double *execution_time, char *params_type, ...){
	va_list variadic_list;
	cl_context _context;
	cl_command_queue _queue;
	cl_mem _memory;
	cl_kernel _kernel;
	cl_program _program;
	cl_event *_profiling_event = (cl_event*)malloc(sizeof(cl_event) * runs);
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[1], local_work_size[1];
	
	//Pre analysis
	param_element = strlen(params_type);	
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	if(verbose){
		printf("Run kernel function: %s\nDevice:\n", function);
		fbclPrintDeviceInfo(device);
	}
	
	_context = clCreateContext(NULL, 1, &device, NULL, NULL, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("\nDEBUG: context\t%s\n", message);
	}
	
	
	printf("Creating queue\n");
	// ONLY ONE QUEUE will be created
	// CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
	// CL_QUEUE_PROFILING_ENABLE
	_queue = clCreateCommandQueue(_context, device, CL_QUEUE_PROFILING_ENABLE, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: queue\t%s\n", message);		
	}
	
	_program = clCreateProgramWithSource(_context, 1, (const char **)&source, (const size_t *)&source_size, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: program\t%s\n", message);
	}

	//Build
	_ret = clBuildProgram(_program, 1, &device, NULL, NULL, NULL);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: build\t%s\n", message);
		fbclProgramBuildInfo(_program, device);
	}
	
	_kernel = clCreateKernel(_program, function, &_ret);
	if(verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: kernel\t%s\n", message);
	}
		
	//Memory INIT
	for (i = 0; i < param_element; i++){
		params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{
		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		if(verbose){ printf("Param %i: Type %c Size %ld\n", i, params_type[i], params[i].param_size); }
		params[i].memory = clCreateBuffer(_context, params[i].memory_type, params[i].param_size, NULL, &_ret);
		if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: memory\t%s\n", message); }
	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message); }
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		_ret = clSetKernelArg(_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: init\t%s\n", message); }
	}


	//Execute
	
	global_work_size[0] = work_items;
	local_work_size[0] = work_items;
	
	for(i=0; i<runs; i++){
		_ret = clEnqueueNDRangeKernel(_queue, _kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event[i]);
		
		if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message); }
		clFlush(_queue); //?
	}
	clFinish(_queue);		// not needed, since we will wait for the event termination
												// see http://stackoverflow.com/questions/11363780/how-to-profile-sequential-launched-multiple-opencl-kernels-by-one-clfinish
	
	printf("Queue finished\n");
	
	cl_int ret;
	
	/*
	if((ret = clWaitForEvents(1, &_profiling_event[0])) != CL_SUCCESS){
		printf("Something wrong with the clWaitForEvent function. Error: %d\n", ret);
	}
	if((ret = clWaitForEvents(1, &_profiling_event[runs-1])) != CL_SUCCESS){
		printf("Something wrong with the clWaitForEvent function. Error: %d\n", ret);
	}
	 */
	
	
	
	/*
	for(i=0; i<runs; i++){
		if((ret = clWaitForEvents(1, &_profiling_event[i])) != CL_SUCCESS){
			printf("Something wrong with the clWaitForEvent function. Error: %d\n", ret);
		}
	}*/
	
	
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event[0],CL_PROFILING_COMMAND_START,sizeof(cl_ulong),&start,NULL);
	clGetEventProfilingInfo(_profiling_event[runs-1],CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end,NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);
	
	
	
	
	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if(verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message); }
		}
	}

	//Destroy all
	_ret = clFlush(_queue);
	_ret = clFinish(_queue);
	for(i=0; i<runs; i++) _ret = clReleaseEvent(_profiling_event[i]);
	_ret = clReleaseKernel(_kernel);
	_ret = clReleaseProgram(_program);	
	for (i = 0; i < param_element; i++){
		_ret = clReleaseMemObject(params[i].memory);
	}
	_ret = clReleaseCommandQueue(_queue);
	_ret = clReleaseContext(_context);
	free(params);

	if(verbose) printf("DEBUG: finish\n");
	return;
}



void fbclInitKernel(cl_device_id device, cl_context *_context, cl_command_queue *_queue, char * source, size_t source_size, cl_program *_program, cl_kernel *_kernel, char * function, int verbose){
	char message[MEM_SIZE];
	cl_int _ret;

	//cl_cont

	*_context = clCreateContext(NULL, 1, &device, NULL, NULL, &_ret);
	if (verbose){
		fbclErrorTranslate(_ret, message);
		printf("\nDEBUG: context\t%s\n", message);
	}

	*_queue = clCreateCommandQueue(*_context, device, CL_QUEUE_PROFILING_ENABLE, &_ret);
	if (verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: queue\t%s\n", message);
	}

	*_program = clCreateProgramWithSource(*_context, 1, (const char **)&source, (const size_t *)&source_size, &_ret);
	if (verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: program\t%s\n", message);
	}

	//Build
	_ret = clBuildProgram(*_program, 1, &device, NULL, NULL, NULL);
	if (verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: build\t%s\n", message);
		fbclProgramBuildInfo(*_program, device);
	}

	*_kernel = clCreateKernel(*_program, function, &_ret);
	if (verbose){
		fbclErrorTranslate(_ret, message);
		printf("DEBUG: kernel\t%s\n", message);
	}
}

int fbclRunSequential(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_program *_program, cl_kernel *_kernel, int work_items, char *source, size_t source_size, char *function, int verbose, cl_double *execution_time, char *params_type, ...){
	va_list variadic_list;
	cl_mem _memory;
	cl_event _profiling_event;
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[1], local_work_size[1];
	int local_mem_used = 0;


	//Pre analysis
	param_element = strlen(params_type);
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	if (verbose){
		printf("Run kernel function: %s\nDevice:\n", function);
		fbclPrintDeviceInfo(_device);
	}

	

	//Memory INIT
	for (i = 0; i < param_element; i++){
		params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{

		// NB: R/W are intended by the device (read-only: DEVICE will only read from it)

		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		case 'l':
			params[i].memory_type = CL_MEM_TYPE_LOCAL;		// custom definition: see fbcl.h
			local_mem_used += params[i].param_size;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		if (verbose){
			printf("Param %i: Type %c Size %ld\n", i, params_type[i], params[i].param_size);
		}

		
		
		if(params[i].memory_type != CL_MEM_TYPE_LOCAL){
			params[i].memory = clCreateBuffer(*_context, params[i].memory_type, params[i].param_size, NULL, &_ret);

			if (verbose){
				fbclErrorTranslate(_ret, message);
				printf("DEBUG: memory\t%s\n", message);
			}
		}
	}

	if(verbose) printf("\n** INFO: kernel will use %d bytes of local memory.\n\n", local_mem_used);
	cl_ulong device_local_memory = 0;
	clGetDeviceInfo(_device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_local_memory, NULL);

	if(local_mem_used > device_local_memory){
		printf("\n\n** ERROR: kernel wants to allocate more than maximum local memory (%ld).\n", device_local_memory);
		printf("\tSome of the parameters may not be allocated as local.\n\n");

		// Select which parameter to dismiss as local
		// TO DO
		/*
		_ret = clReleaseKernel(*_kernel);
		_ret = clReleaseProgram(*_program);
		for (i = 0; i < param_element; i++){
			if (params[i].memory_type != CL_MEM_TYPE_LOCAL)
				_ret = clReleaseMemObject(params[i].memory);
		}
		_ret = clReleaseCommandQueue(*_queue);
		_ret = clReleaseContext(*_context);
		*/
		free(params);

		printf("Exiting...\n");
		return -1;

	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message); }
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		if(params[i].memory_type == CL_MEM_TYPE_LOCAL){
			if(verbose) printf("DEBUG: setting local parameter #%d with size %ld\n", i, params[i].param_size);
			_ret = clSetKernelArg(*_kernel, i, (size_t) params[i].param_size, NULL);
		}else{
			_ret = clSetKernelArg(*_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		}
		if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: init\t%s\n", message); }
	}


	//Execute
	//_ret = clEnqueueTask(_queue, _kernel, 0, NULL, NULL);
	global_work_size[0] = work_items;
	local_work_size[0] = work_items;
	_ret = clEnqueueNDRangeKernel(*_queue, *_kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event);
	if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message); }

	clWaitForEvents(1, &_profiling_event);
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);

	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message); }
		}
	}

	
	//Destroy only run-specific structures
	
	_ret = clFlush(*_queue);
	_ret = clFinish(*_queue);
	_ret = clReleaseEvent(_profiling_event);
	for (i = 0; i < param_element; i++){
		if(params[i].memory_type != CL_MEM_TYPE_LOCAL){
			_ret = clReleaseMemObject(params[i].memory);
		}
	}
	free(params);

	if (verbose) printf("DEBUG: finish\n");
	return 0;
}


int fbclRunParallelSingle(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_kernel *_kernel,	int work_items, char *function, int verbose, cl_double *execution_time, char *params_type, ...)
{
	va_list variadic_list;
	cl_event _profiling_event;
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[2], local_work_size[2];
	int local_mem_used = 0;


	//Pre analysis
	param_element = strlen(params_type);
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	if (verbose){
		printf("Run kernel function: %s\nDevice:\n", function);
		fbclPrintDeviceInfo(_device);
	}



	//Memory INIT
	for (i = 0; i < param_element; i++){
		params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{

			// NB: R/W are intended by the device (read-only: DEVICE will only read from it)

		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		case 'l':
			params[i].memory_type = CL_MEM_TYPE_LOCAL;		// custom definition: see fbcl.h
			local_mem_used += params[i].param_size;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		if (verbose){
			printf("Param %i: Type %c Size %ld\n", i, params_type[i], params[i].param_size);
		}



		if (params[i].memory_type != CL_MEM_TYPE_LOCAL){
			params[i].memory = clCreateBuffer(*_context, params[i].memory_type, params[i].param_size, NULL, &_ret);

			if (verbose){
				fbclErrorTranslate(_ret, message);
				printf("DEBUG: memory\t%s\n", message);
			}
		}
	}

	if (verbose) printf("\n** INFO: kernel will use %d bytes of local memory.\n\n", local_mem_used);
	cl_ulong device_local_memory = 0;
	clGetDeviceInfo(_device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_local_memory, NULL);

	if (local_mem_used > device_local_memory){
		printf("\n\n** ERROR: kernel wants to allocate more than maximum local memory (%ld).\n", device_local_memory);
		printf("\tSome of the parameters may not be allocated as local.\n\n");

		free(params);

		printf("Exiting...\n");
		return -1;

	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message); }
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_TYPE_LOCAL){
			if (verbose) printf("DEBUG: setting local parameter #%d with size %ld\n", i, params[i].param_size);
			_ret = clSetKernelArg(*_kernel, i, (size_t)params[i].param_size, NULL);
		}
		else{
			_ret = clSetKernelArg(*_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		}
		if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: init\t%s\n", message); }
	}


	/*********************************************************   Execute   *************************************************************/
	global_work_size[0] = LINEAR_DIMENSION;
	//global_work_size[1] = VERTICAL_GROUPS;
	global_work_size[1] = 0;
	
	local_work_size[0] = work_items;
	//local_work_size[1] = 1;
	local_work_size[1] = 0;
	printf("RUN KERNEL with global_size[0]=%ld and local_size[0]=%ld => horizontal group number=%ld/%ld=%ld\n", global_work_size[0], local_work_size[0], global_work_size[0], local_work_size[0], global_work_size[0]/local_work_size[0]);
	printf("Vertical dimensions = global_work_size[1]=%ld, local_work_size[1]=%ld\n\n", global_work_size[1], local_work_size[1]);
	
	_ret = clEnqueueNDRangeKernel(*_queue, *_kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event);
	
	if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message); }

	clWaitForEvents(1, &_profiling_event);
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);

	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message); }
		}
	}


	//Destroy only run-specific structures

	_ret = clFlush(*_queue);
	_ret = clFinish(*_queue);
	_ret = clReleaseEvent(_profiling_event);
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type != CL_MEM_TYPE_LOCAL){
			_ret = clReleaseMemObject(params[i].memory);
		}
	}
	free(params);

	if (verbose) printf("DEBUG: finish\n");
	return 0;
}

int fbclRunParallel(cl_device_id _device, cl_context *_context, cl_command_queue *_queue, cl_kernel *_kernel, int work_groups, int work_items, char *function, int verbose, cl_double *execution_time, char *params_type, ...){
	va_list variadic_list;
	cl_event _profiling_event;
	cl_int _ret;
	int i, param_element;
	FBCL_PARAM *params;
	char message[MEM_SIZE];
	size_t global_work_size[2], local_work_size[2];
	int local_mem_used = 0;


	//Pre analysis
	param_element = strlen(params_type);
	va_start(variadic_list, params_type);
	params = (FBCL_PARAM *)malloc(sizeof(FBCL_PARAM)*param_element);

	//INIT
	if (verbose){
		printf("Run kernel function: %s\nDevice:\n", function);
		fbclPrintDeviceInfo(_device);
	}



	//Memory INIT
	for (i = 0; i < param_element; i++){
		//params[i].param_size = (int)va_arg(variadic_list, int);
		params[i].param_size = (unsigned long)va_arg(variadic_list, unsigned long);
		params[i].param = (void*)va_arg(variadic_list, void*);
		switch (params_type[i])
		{

			// NB: R/W are intended by the device (read-only: DEVICE will only read from it)

		case 'a':
			params[i].memory_type = CL_MEM_READ_WRITE;
			break;
		case 'w':
			params[i].memory_type = CL_MEM_WRITE_ONLY;
			break;
		case 'r':
			params[i].memory_type = CL_MEM_READ_ONLY;
			break;
		case 'l':
			params[i].memory_type = CL_MEM_TYPE_LOCAL;		// custom definition: see fbcl.h
			local_mem_used += params[i].param_size;
			break;
		default:
			params[i].memory_type = CL_MEM_READ_WRITE;
		}
		if (verbose){
			printf("Param %i: Type %c Size %ld\n", i, params_type[i], params[i].param_size);
		}



		if (params[i].memory_type != CL_MEM_TYPE_LOCAL){
			params[i].memory = clCreateBuffer(*_context, params[i].memory_type, (size_t)params[i].param_size, NULL, &_ret);

			if (verbose){
				fbclErrorTranslate(_ret, message);
				printf("DEBUG: memory\t%s\n", message);
			}
		}
	}

	if (verbose) printf("\n** INFO: kernel will use %d bytes of local memory.\n\n", local_mem_used);
	cl_ulong device_local_memory = 0;
	clGetDeviceInfo(_device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_local_memory, NULL);

	if (local_mem_used > device_local_memory){
		printf("\n\n** ERROR: kernel wants to allocate more than maximum local memory (%ld).\n", device_local_memory);
		printf("\tSome of the parameters may not be allocated as local.\n\n");

		free(params);

		printf("Exiting...\n");
		return -1;

	}

	//Put data for read and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_READ_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueWriteBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			//_ret = clEnqueueReadBuffer(_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: write\t%s\n", message); }
		}
	}

	//Kernel INIT
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_TYPE_LOCAL){
			if (verbose) printf("DEBUG: setting local parameter #%d with size %ld\n", i, params[i].param_size);
			_ret = clSetKernelArg(*_kernel, i, (size_t)params[i].param_size, NULL);
		}
		else{
			_ret = clSetKernelArg(*_kernel, i, sizeof(cl_mem), (void *)&(params[i].memory));
		}
		if (verbose){
			fbclErrorTranslate(_ret, message);
			printf("DEBUG: init\t%s\n", message);
			if(strcmp(message, "CL_SUCCESS") != 0)
				printf("\tparam_size = %ld\n", params[i].param_size);
		}
	}


	/*********************************************************   Execute   *************************************************************/
	
	local_work_size[0] = work_items;
	local_work_size[1] = 1;
	
	
	global_work_size[0] = MIN(LINEAR_DIMENSION, work_groups*work_items);
	
	
	int wg_per_lin = global_work_size[0] / local_work_size[0];
	global_work_size[1] = MIN(LINEAR_DIMENSION, work_groups / wg_per_lin + 1);
	
	
	printf("\n#RUN KERNEL with global_size[0]=%ld and local_size[0]=%ld => horizontal group number=%ld/%ld=%ld\n", global_work_size[0], local_work_size[0], global_work_size[0], local_work_size[0], global_work_size[0]/local_work_size[0]);
	printf("Vertical dimensions = global_work_size[1]=%ld, local_work_size[1]=%ld\n\n", global_work_size[1], local_work_size[1]);
	printf("\n\n * parallel=%d, workers=%d ==> hor_groups=%ld, vertical_dim=%ld\n\n\n\n\n", work_groups, work_items, global_work_size[0]/local_work_size[0], global_work_size[1]);
	
		
	_ret = clEnqueueNDRangeKernel(*_queue, *_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, &_profiling_event);
	
	if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: run\t%s\n", message); }

	clWaitForEvents(1, &_profiling_event);
	cl_ulong start = 0, end = 0;
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(_profiling_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	*execution_time = (double)(end - start)*(double)(1e-06);

	//Get data for write and read-write
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type == CL_MEM_WRITE_ONLY || params[i].memory_type == CL_MEM_READ_WRITE){
			_ret = clEnqueueReadBuffer(*_queue, params[i].memory, CL_TRUE, 0, params[i].param_size, params[i].param, 0, NULL, NULL);
			if (verbose){ fbclErrorTranslate(_ret, message); printf("DEBUG: read\t%s\n", message); }
		}
	}


	//Destroy only run-specific structures

	_ret = clFlush(*_queue);
	_ret = clFinish(*_queue);
	_ret = clReleaseEvent(_profiling_event);
	for (i = 0; i < param_element; i++){
		if (params[i].memory_type != CL_MEM_TYPE_LOCAL){
			_ret = clReleaseMemObject(params[i].memory);
		}
	}
	free(params);

	if (verbose) printf("DEBUG: finish\n");
	return 0;
}

int fbclDestroy(cl_context *_context, cl_program *_program, cl_kernel *_kernel, cl_command_queue *_queue){
	int _ret;
	_ret = clReleaseKernel(*_kernel);
	_ret = clReleaseProgram(*_program);
	_ret = clReleaseCommandQueue(*_queue);
	_ret = clReleaseContext(*_context);

	return _ret;
}


char *fbclGetSourceFromFile(const char *file_name, size_t *source_size){
	FILE *fp;
	long long fp_size;
	char *source_str;
	
	fp = fopen(file_name,"r");
	if (!fp) {
		printf("ERR: Failed to load kernel file\n");
		exit(1);
	}
	fseek(fp, 0, SEEK_END);
	fgetpos(fp, (fpos_t*)&fp_size);
	fseek(fp, 0, SEEK_SET);

	source_str = (char*)malloc((long long)sizeof(char)*fp_size);
	*source_size = fread(source_str, 1, fp_size, fp);
	if (ferror(fp) != 0){
		printf("ERR: Failed to read kernel file\n");
		exit(1);
	}
	fclose(fp);
	return source_str;
}

int fbclGetFirstPlatform(cl_platform_id *platform){
	cl_platform_id platforms[1];
	if (fbclNumberOfPlatforms() == 0)
		return 0;
	clGetPlatformIDs(1, platforms, NULL);
	*platform = platforms[0];
	return 1;
}

int fbclGetFirstDevice(cl_device_id *device, cl_device_type type){
	cl_platform_id platform;
	cl_device_id devices[1];
	if (!fbclGetFirstPlatform(&platform))
		return 0;
	if (fbclNumberOfDevices(platform, type) == 0)
		return 0;
	clGetDeviceIDs(platform, type, 1, devices, NULL);
	*device = devices[0];
	return 1;
}

int fbclGetFirstGPU(cl_device_id *device){
	return fbclGetFirstDevice(device, CL_DEVICE_TYPE_GPU);
}

int fbclGetFirstCPU(cl_device_id *device){
	return fbclGetFirstDevice(device, CL_DEVICE_TYPE_CPU);
}

cl_uint fbclNumberOfPlatforms(){
	cl_uint n;
	clGetPlatformIDs(0, NULL, &n);
	return n;
}

cl_uint fbclNumberOfDevices(cl_platform_id platform, cl_device_type type){
	cl_uint n;
	clGetDeviceIDs(platform, type, 0, NULL, &n);
	return n;
}

void fbclDebug(){
	int ret, queryInt;
	char buffer[1024];
	cl_platform_id *platforms;
	cl_device_id *devices;
	cl_uint num_platforms;
	cl_uint num_devices;
	clGetPlatformIDs(0, NULL, &num_platforms);
	platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id)*num_platforms);
	clGetPlatformIDs(num_platforms, platforms, NULL);
	for (int i = 0; i < num_platforms; i++){
		clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 1024, buffer, NULL);
		printf("Platform %d: %s\n", i, buffer);

		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
		devices = (cl_device_id *)malloc(sizeof(cl_device_id)*num_devices);
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);

		for (int j = 0; j < num_devices; j++){
			printf("\n");
			fbclPrintDeviceInfo(devices[j]);
			/*
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 1024, buffer, NULL);
			printf("\tDevice %d: %s\n", j, buffer);
			clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, 1024, buffer, NULL);
			printf("\t\tVendor: %s\n", buffer);
			clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 1024, buffer, NULL);
			printf("\t\tDeviceVersion %s\n", buffer);
			clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 1024, buffer, NULL);
			printf("\t\tDriverVersion %s\n", buffer);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(int), &queryInt, NULL);
			printf("\t\tComputeUnit %d\n", queryInt); 
			*/
		}

	}
	return;
}

void fbclGetDeviceLocalMemory(cl_device_id device, cl_ulong * local_memory_size, cl_int * error){
	*error = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), local_memory_size, NULL);
	return;
}

void fbclGetMaxWorkGroupSize(cl_device_id device, size_t *size, cl_int *error, char * message){
	char * buf = (char*) calloc(4096, sizeof(char));
	*error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), size, NULL);
	fbclErrorTranslate(*error, buf);
	
	if(message != NULL) strcpy(message, buf);
	free(buf);
	return;
}

void fbclPrintDeviceInfo(cl_device_id device){
	char queryBuffer[1024];
	size_t sizes[3];
	cl_int clError;
	cl_device_local_mem_type lmtype;
	cl_uint ret;
	cl_ulong l;
	
	clError = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("\tDeviceName: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	clError = clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("\tDeviceVendor: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	clError = clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("\tDeviceVersion: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	clError = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("\tDriverVersion: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	clError = clGetDeviceInfo(device, CL_DEVICE_VENDOR_ID, sizeof(cl_uint), &ret, NULL);
	printf("\tVendor ID: %d\n", ret);
	clError = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &ret, NULL);
	printf("\tMax Compute Units: %d\n", ret);
	clError = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &ret, NULL);
	printf("\tMax work-item dimensions: %d\n", ret);
	clError = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, 3*sizeof(size_t), sizes, NULL);
	printf("\tMax work-item sizes: (%ld, %ld, %ld) \n", sizes[0], sizes[1], sizes[2]);
	clError = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &(sizes[0]), NULL);
	printf("\tMax work-group size: %ld\n", sizes[0]);
	clError = clGetDeviceInfo(device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &(sizes[0]), NULL);
	printf("\tPreferred work-group size: %ld\n", sizes[0]);
	clError = clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(size_t), &(sizes[0]), NULL);
	printf("\tImage support: %ld ", sizes[0]); if(sizes[0]){ printf("(Yes)\n"); }else{ printf("(No)\n"); }
	clError = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &l, NULL);
	printf("\tDevice Local Memory size: %ld\n", l);
	clError = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &lmtype, NULL);
	printf("\tDevice Local Memory type: ");
		if(lmtype == CL_LOCAL)
			printf("CL_LOCAL\n");
		else if(lmtype == CL_GLOBAL)
			printf("CL_GLOBAL\n");
		else
			printf("Unknown type\n");

}


void fbclErrorTranslate(cl_int value, char *message){	
	char values[][50] = { 
		"CL_SUCCESS",
		"CL_DEVICE_NOT_FOUND",
		"CL_DEVICE_NOT_AVAILABLE",
		"CL_COMPILER_NOT_AVAILABLE",
		"CL_MEM_OBJECT_ALLOCATION_FAILURE",
		"CL_OUT_OF_RESOURCES",
		"CL_OUT_OF_HOST_MEMORY",
		"CL_PROFILING_INFO_NOT_AVAILABLE",
		"CL_MEM_COPY_OVERLAP",
		"CL_IMAGE_FORMAT_MISMATCH",
		"CL_IMAGE_FORMAT_NOT_SUPPORTED",
		"CL_BUILD_PROGRAM_FAILURE",
		"CL_MAP_FAILURE",
		"CL_INVALID_VALUE",
		"CL_INVALID_DEVICE_TYPE",
		"CL_INVALID_PLATFORM",
		"CL_INVALID_DEVICE",
		"CL_INVALID_CONTEXT",
		"CL_INVALID_QUEUE_PROPERTIES",
		"CL_INVALID_COMMAND_QUEUE",
		"CL_INVALID_HOST_PTR",
		"CL_INVALID_MEM_OBJECT",
		"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
		"CL_INVALID_IMAGE_SIZE",
		"CL_INVALID_SAMPLER",
		"CL_INVALID_BINARY",
		"CL_INVALID_BUILD_OPTIONS",
		"CL_INVALID_PROGRAM",
		"CL_INVALID_PROGRAM_EXECUTABLE",
		"CL_INVALID_KERNEL_NAME",
		"CL_INVALID_KERNEL_DEFINITION",
		"CL_INVALID_KERNEL",
		"CL_INVALID_ARG_INDEX",
		"CL_INVALID_ARG_VALUE",
		"CL_INVALID_ARG_SIZE",
		"CL_INVALID_KERNEL_ARGS",
		"CL_INVALID_WORK_DIMENSION",
		"CL_INVALID_WORK_GROUP_SIZE",
		"CL_INVALID_WORK_ITEM_SIZE",
		"CL_INVALID_GLOBAL_OFFSET",
		"CL_INVALID_EVENT_WAIT_LIST",
		"CL_INVALID_EVENT",
		"CL_INVALID_OPERATION",
		"CL_INVALID_GL_OBJECT",
		"CL_INVALID_BUFFER_SIZE",
		"CL_INVALID_MIP_LEVEL",
		"CL_INVALID_GLOBAL_WORK_SIZE"};
	value = abs(value);
	if (value > 12)
		value -= 17;
	strcpy(message, values[value]);
	return;
}

void fbclProgramBuildInfo(cl_program program, cl_device_id device){
	cl_build_status status;
	char *log_str;
	size_t log_size;
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
	if (status == CL_BUILD_NONE){
		printf("BUILD: No Info\n");
	}
	else if (status == CL_BUILD_ERROR){
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		log_str = (char*)malloc(sizeof(char)*log_size);
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(char)*log_size, log_str, NULL);
		printf("BUILD: %s\n", log_str);
		free(log_str);
	}
	return;
}

void fbclVariadicGeneral(long int *somma, int *massimo, int *minimo, int count, ...){
	va_list list;
	long int _sum = 0;
	int value, _max, _min;
	va_start(list, count);
	_max = INT_MIN;
	_min = INT_MAX;
	for (int i = 0; i<count; i++){
		value = va_arg(list, int);
		_sum += value;
		if (value>_max)
			_max = value;
		if (value<_min)
			_min = value;
	}
	
	if (!(somma == NULL))
		*somma = _sum;
	
	if (!(minimo == NULL))
		*minimo = _min;
	
	if (!(massimo == NULL))
		*massimo = _max;
	return;
}

void fbclVariadicGeneral16(long int *somma, int *massimo, int *minimo, int count, ...){
	va_list list;
	long int _sum = 0;
	int value, _max, _min;
	va_start(list, count);
	_max = INT_MIN;
	_min = INT_MAX;
	for (int i = 0; i<count; i++){
		value = va_arg(list, int);
		_sum += value;
		if (value>_max)
			_max = value;
		if (value<_min)
			_min = value;
	}


	if (!(somma == NULL))
		*somma = _sum;

	if (_min % 16){
		_min += 16 - _min % 16;
	}
	if (!(minimo == NULL))
		*minimo = _min;

	if (_max % 16){
		_max += 16 - _max % 16;
	}
	if (!(massimo == NULL))
		*massimo = _max;
	return;
}


void fbclVariadicGeneral32(long int *somma, int *massimo, int *minimo, int count, ...){
	va_list list;
	long int _sum = 0;
	int value, _max, _min;
	va_start(list, count);
	_max = INT_MIN;
	_min = INT_MAX;
	for (int i = 0; i<count; i++){
		value = va_arg(list, int);
		_sum += value;
		if (value>_max)
			_max = value;
		if (value<_min)
			_min = value;
	}


	if (!(somma == NULL))
		*somma = _sum;

	if (_min % 32){
		_min += 32 - _min % 32;
	}
	if (!(minimo == NULL))
		*minimo = _min;

	if (_max % 32){
		_max += 32 - _max % 32;
	}
	if (!(massimo == NULL))
		*massimo = _max;
	return;
}
