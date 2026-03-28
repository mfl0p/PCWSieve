/* #######################################################################
    Copyright 2011 Oscar Amoros Huguet, Cristian Garcia Marin

    This file is part of SimpleOpenCL

    SimpleOpenCL is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    SimpleOpenCL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SimpleOpenCL. If not, see <http://www.gnu.org/licenses/>.

   ####################################################################### 

   SimpleOpenCL Version 0.010_27_02_2013 

*/

#ifdef __cplusplus
extern "C" {
#endif

#include "simpleCL.h"

void sclPrintErrorFlags( cl_int flag ){
    
	switch (flag){

		case CL_DEVICE_NOT_FOUND:
			printf("\nCL_DEVICE_NOT_FOUND\n");
			break;
		case CL_DEVICE_NOT_AVAILABLE:
			printf("\nCL_DEVICE_NOT_AVAILABLE\n");
			break;
		case CL_COMPILER_NOT_AVAILABLE:
			printf("\nCL_COMPILER_NOT_AVAILABLE\n");
			break;
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			printf("\nCL_PROFILING_INFO_NOT_AVAILABLE\n");
			break;
		case CL_MEM_COPY_OVERLAP:
			printf("\nCL_MEM_COPY_OVERLAP\n");
			break;
		case CL_IMAGE_FORMAT_MISMATCH:
			printf("\nCL_IMAGE_FORMAT_MISMATCH\n");
			break;
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			printf("\nCL_IMAGE_FORMAT_NOT_SUPPORTED\n");
			break;
		case CL_INVALID_COMMAND_QUEUE:
			printf("\nCL_INVALID_COMMAND_QUEUE\n");
			break;
		case CL_INVALID_CONTEXT:
			printf("\nCL_INVALID_CONTEXT\n");
			break;
		case CL_INVALID_MEM_OBJECT:
			printf("\nCL_INVALID_MEM_OBJECT\n");
			break;
		case CL_INVALID_VALUE:
			printf("\nCL_INVALID_VALUE\n");
			break;
		case CL_INVALID_EVENT_WAIT_LIST:
			printf("\nCL_INVALID_EVENT_WAIT_LIST\n");
			break;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			printf("\nCL_MEM_OBJECT_ALLOCATION_FAILURE\n");
			break;
		case CL_OUT_OF_HOST_MEMORY:
			printf("\nCL_OUT_OF_HOST_MEMORY\n");
			break;

		case CL_INVALID_PROGRAM_EXECUTABLE:
			printf("\nCL_INVALID_PROGRAM_EXECUTABLE\n");
			break;
		case CL_INVALID_KERNEL:
			printf("\nCL_INVALID_KERNEL\n");
			break;
		case CL_INVALID_KERNEL_ARGS:
			printf("\nCL_INVALID_KERNEL_ARGS\n");
			break;
		case CL_INVALID_WORK_DIMENSION:
			printf("\nCL_INVALID_WORK_DIMENSION\n");
			break;
#ifndef __APPLE__ 
		case CL_INVALID_GLOBAL_WORK_SIZE:
			printf("\nCL_INVALID_GLOBAL_WORK_SIZE\n");
			break;
#endif
		case CL_INVALID_WORK_GROUP_SIZE:
			printf("\nCL_INVALID_WORK_GROUP_SIZE\n");
			break;
		case CL_INVALID_WORK_ITEM_SIZE:
			printf("\nCL_INVALID_WORK_ITEM_SIZE\n");
			break;
		case CL_INVALID_GLOBAL_OFFSET:
			printf("\nCL_INVALID_GLOBAL_OFFSET\n");
			break;
		case CL_OUT_OF_RESOURCES:
			printf("\nCL_OUT_OF_RESOURCES\n");
			break;

		case CL_INVALID_PROGRAM:
			printf("\nCL_INVALID_PROGRAM\n");
			break;
		case CL_INVALID_KERNEL_NAME:
			printf("\nCL_INVALID_KERNEL_NAME\n");
			break;
		case CL_INVALID_KERNEL_DEFINITION:
			printf("\nCL_INVALID_KERNEL_DEFINITION\n");
			break;
		case CL_INVALID_BUFFER_SIZE:
			printf("\nCL_INVALID_BUFFER_SIZE\n");
			break;
		case CL_BUILD_PROGRAM_FAILURE:
			printf("\nCL_BUILD_PROGRAM_FAILURE\n");
			break;
		case CL_INVALID_ARG_INDEX:
			printf("\nCL_INVALID_ARG_INDEX\n");
			break;
		case CL_INVALID_ARG_VALUE:
			printf("\nCL_INVALID_ARG_VALUE\n");
			break;
		case CL_MAP_FAILURE:
			printf("\nCL_MAP_FAILURE\n");
			break;
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			printf("\nCL_MISALIGNED_SUB_BUFFER_OFFSET\n");
			break;
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			printf("\nCL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST\n");
			break;
		case CL_INVALID_DEVICE_TYPE:
			printf("\nCL_INVALID_DEVICE_TYPE\n");
			break;
		case CL_INVALID_PLATFORM:
			printf("\nCL_INVALID_PLATFORM\n");
			break;
		case CL_INVALID_DEVICE:
			printf("\nCL_INVALID_DEVICE\n");
			break; 
		case CL_INVALID_QUEUE_PROPERTIES:
			printf("\nCL_INVALID_QUEUE_PROPERTIES\n");
			break; 
		case CL_INVALID_HOST_PTR:
			printf("\nCL_INVALID_HOST_PTR\n");
			break;
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			printf("\nCL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
			break;
		case CL_INVALID_IMAGE_SIZE:
			printf("\nCL_INVALID_IMAGE_SIZE\n");
			break;
		case CL_INVALID_SAMPLER:
			printf("\nCL_INVALID_SAMPLER\n");
			break;
		case CL_INVALID_BINARY:
			printf("\nCL_INVALID_BINARY\n");
			break;
		case CL_INVALID_BUILD_OPTIONS:
			printf("\nCL_INVALID_BUILD_OPTIONS\n");
			break;
		case CL_INVALID_ARG_SIZE:
			printf("\nCL_INVALID_ARG_SIZE\n");
			break;
		case CL_INVALID_EVENT:
			printf("\nCL_INVALID_EVENT\n");
			break;
		case CL_INVALID_OPERATION:
			printf("\nCL_INVALID_OPERATION\n");
			break;
		case CL_INVALID_GL_OBJECT:
			printf("\nCL_INVALID_GL_OBJECT\n");
			break;
		case CL_INVALID_MIP_LEVEL:
			printf("\nCL_INVALID_MIP_LEVEL\n");
			break;
		case CL_INVALID_PROPERTY:
			printf("\nCL_INVALID_PROPERTY\n");
			break;
		default:
			printf("\nUnknown error code: %d\n",flag);    
	}


	switch (flag){

		case CL_DEVICE_NOT_FOUND:
			fprintf(stderr, "\nCL_DEVICE_NOT_FOUND\n");
			break;
		case CL_DEVICE_NOT_AVAILABLE:
			fprintf(stderr, "\nCL_DEVICE_NOT_AVAILABLE\n");
			break;
		case CL_COMPILER_NOT_AVAILABLE:
			fprintf(stderr, "\nCL_COMPILER_NOT_AVAILABLE\n");
			break;
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			fprintf(stderr, "\nCL_PROFILING_INFO_NOT_AVAILABLE\n");
			break;
		case CL_MEM_COPY_OVERLAP:
			fprintf(stderr, "\nCL_MEM_COPY_OVERLAP\n");
			break;
		case CL_IMAGE_FORMAT_MISMATCH:
			fprintf(stderr, "\nCL_IMAGE_FORMAT_MISMATCH\n");
			break;
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			fprintf(stderr, "\nCL_IMAGE_FORMAT_NOT_SUPPORTED\n");
			break;
		case CL_INVALID_COMMAND_QUEUE:
			fprintf(stderr, "\nCL_INVALID_COMMAND_QUEUE\n");
			break;
		case CL_INVALID_CONTEXT:
			fprintf(stderr, "\nCL_INVALID_CONTEXT\n");
			break;
		case CL_INVALID_MEM_OBJECT:
			fprintf(stderr, "\nCL_INVALID_MEM_OBJECT\n");
			break;
		case CL_INVALID_VALUE:
			fprintf(stderr, "\nCL_INVALID_VALUE\n");
			break;
		case CL_INVALID_EVENT_WAIT_LIST:
			fprintf(stderr, "\nCL_INVALID_EVENT_WAIT_LIST\n");
			break;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			fprintf(stderr, "\nCL_MEM_OBJECT_ALLOCATION_FAILURE\n");
			break;
		case CL_OUT_OF_HOST_MEMORY:
			fprintf(stderr, "\nCL_OUT_OF_HOST_MEMORY\n");
			break;

		case CL_INVALID_PROGRAM_EXECUTABLE:
			fprintf(stderr, "\nCL_INVALID_PROGRAM_EXECUTABLE\n");
			break;
		case CL_INVALID_KERNEL:
			fprintf(stderr, "\nCL_INVALID_KERNEL\n");
			break;
		case CL_INVALID_KERNEL_ARGS:
			fprintf(stderr, "\nCL_INVALID_KERNEL_ARGS\n");
			break;
		case CL_INVALID_WORK_DIMENSION:
			fprintf(stderr, "\nCL_INVALID_WORK_DIMENSION\n");
			break;
#ifndef __APPLE__ 
		case CL_INVALID_GLOBAL_WORK_SIZE:
			fprintf(stderr, "\nCL_INVALID_GLOBAL_WORK_SIZE\n");
			break;
#endif
		case CL_INVALID_WORK_GROUP_SIZE:
			fprintf(stderr, "\nCL_INVALID_WORK_GROUP_SIZE\n");
			break;
		case CL_INVALID_WORK_ITEM_SIZE:
			fprintf(stderr, "\nCL_INVALID_WORK_ITEM_SIZE\n");
			break;
		case CL_INVALID_GLOBAL_OFFSET:
			fprintf(stderr, "\nCL_INVALID_GLOBAL_OFFSET\n");
			break;
		case CL_OUT_OF_RESOURCES:
			fprintf(stderr, "\nCL_OUT_OF_RESOURCES\n");
			break;

		case CL_INVALID_PROGRAM:
			fprintf(stderr, "\nCL_INVALID_PROGRAM\n");
			break;
		case CL_INVALID_KERNEL_NAME:
			fprintf(stderr, "\nCL_INVALID_KERNEL_NAME\n");
			break;
		case CL_INVALID_KERNEL_DEFINITION:
			fprintf(stderr, "\nCL_INVALID_KERNEL_DEFINITION\n");
			break;
		case CL_INVALID_BUFFER_SIZE:
			fprintf(stderr, "\nCL_INVALID_BUFFER_SIZE\n");
			break;
		case CL_BUILD_PROGRAM_FAILURE:
			fprintf(stderr, "\nCL_BUILD_PROGRAM_FAILURE\n");
			break;
		case CL_INVALID_ARG_INDEX:
			fprintf(stderr, "\nCL_INVALID_ARG_INDEX\n");
			break;
		case CL_INVALID_ARG_VALUE:
			fprintf(stderr, "\nCL_INVALID_ARG_VALUE\n");
			break;
		case CL_MAP_FAILURE:
			fprintf(stderr, "\nCL_MAP_FAILURE\n");
			break;
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			fprintf(stderr, "\nCL_MISALIGNED_SUB_BUFFER_OFFSET\n");
			break;
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			fprintf(stderr, "\nCL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST\n");
			break;
		case CL_INVALID_DEVICE_TYPE:
			fprintf(stderr, "\nCL_INVALID_DEVICE_TYPE\n");
			break;
		case CL_INVALID_PLATFORM:
			fprintf(stderr, "\nCL_INVALID_PLATFORM\n");
			break;
		case CL_INVALID_DEVICE:
			fprintf(stderr, "\nCL_INVALID_DEVICE\n");
			break; 
		case CL_INVALID_QUEUE_PROPERTIES:
			fprintf(stderr, "\nCL_INVALID_QUEUE_PROPERTIES\n");
			break; 
		case CL_INVALID_HOST_PTR:
			fprintf(stderr, "\nCL_INVALID_HOST_PTR\n");
			break;
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			fprintf(stderr, "\nCL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
			break;
		case CL_INVALID_IMAGE_SIZE:
			fprintf(stderr, "\nCL_INVALID_IMAGE_SIZE\n");
			break;
		case CL_INVALID_SAMPLER:
			fprintf(stderr, "\nCL_INVALID_SAMPLER\n");
			break;
		case CL_INVALID_BINARY:
			fprintf(stderr, "\nCL_INVALID_BINARY\n");
			break;
		case CL_INVALID_BUILD_OPTIONS:
			fprintf(stderr, "\nCL_INVALID_BUILD_OPTIONS\n");
			break;
		case CL_INVALID_ARG_SIZE:
			fprintf(stderr, "\nCL_INVALID_ARG_SIZE\n");
			break;
		case CL_INVALID_EVENT:
			fprintf(stderr, "\nCL_INVALID_EVENT\n");
			break;
		case CL_INVALID_OPERATION:
			fprintf(stderr, "\nCL_INVALID_OPERATION\n");
			break;
		case CL_INVALID_GL_OBJECT:
			fprintf(stderr, "\nCL_INVALID_GL_OBJECT\n");
			break;
		case CL_INVALID_MIP_LEVEL:
			fprintf(stderr, "\nCL_INVALID_MIP_LEVEL\n");
			break;
		case CL_INVALID_PROPERTY:
			fprintf(stderr, "\nCL_INVALID_PROPERTY\n");
			break;
		default:
			fprintf(stderr, "\nUnknown error code: %d\n",flag);    
	}

	exit(EXIT_FAILURE);
}



void sclGetBinary( sclSoft software ){

	size_t size;
	cl_int err;

	err = clGetProgramInfo( software.program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &size, NULL );
	if ( err!=CL_SUCCESS ) {
		printf( "Error: clGetProgramInfo\n" );
		fprintf(stderr, "Error: clGetProgramInfo\n" );
		sclPrintErrorFlags( err );
	}

	unsigned char * binary = new unsigned char [ size ];

	err = clGetProgramInfo( software.program, CL_PROGRAM_BINARIES, size, &binary, NULL );
	if ( err!=CL_SUCCESS ) {
		printf( "Error: clGetProgramInfo\n" );
		fprintf(stderr, "Error: clGetProgramInfo\n" );
		sclPrintErrorFlags( err );
	}

	FILE * fpbin = fopen( software.kernelName, "wb" );
	if( fpbin == NULL ){
		printf( "Error: cannot write binary: %s\n", software.kernelName );
		fprintf( stderr, "Error: cannot write binary: %s\n", software.kernelName );
	}
	else
	{
		fwrite( binary, 1, size, fpbin );
		fclose( fpbin );
	}
	delete [ ] binary;

}


char* _sclLoadProgramSource( const char *filename )
{ 
	struct stat statbuf;
	FILE *fh; 
	char *source;
	
	fh = fopen( filename, "r" );
	if ( fh == 0 )
		return 0;
	
	stat( filename, &statbuf );
	source = (char *)malloc( statbuf.st_size + 1 );
	size_t x = fread( source, statbuf.st_size, 1, fh );
	if(x != 1){
		printf( "Error on _sclLoadProgramSource\n" );
		fprintf(stderr, "Error on _sclLoadProgramSource\n" );
	}
	source[ statbuf.st_size ] = '\0';
	
	fclose( fh );
	
	return source; 
}
 
cl_program _sclCreateProgram( const char* program_source, cl_context context )
{
	cl_program program;

	cl_int err;
	
	program = clCreateProgramWithSource( context, 1, &program_source, NULL, &err );
	if ( err!=CL_SUCCESS ) {
		printf( "Error on createProgram\n" );
		fprintf(stderr, "Error on createProgram\n" );
		sclPrintErrorFlags( err );
	}

	
	return program;
}

void _sclBuildProgram( cl_program program, cl_device_id devices, const char* pName, const char * options )
{
	cl_int err;
	char build_c[4096];
	
	err = clBuildProgram( program, 0, NULL, options, NULL, NULL );

	// print nvidia kernel build log
//	err = clBuildProgram( program, 0, NULL, "-cl-nv-verbose", NULL, NULL );
//	clGetProgramBuildInfo( program, devices, CL_PROGRAM_BUILD_LOG, 4096, build_c, NULL );
//	printf( "Build Log for %s_program:\n%s\n", pName, build_c );

   	if ( err != CL_SUCCESS ) {
		printf( "Error on buildProgram " );
		fprintf( stdout, "\nRequestingInfo\n" );
		clGetProgramBuildInfo( program, devices, CL_PROGRAM_BUILD_LOG, 4096, build_c, NULL );
		printf( "Build Log for %s_program:\n%s\n", pName, build_c );

		fprintf(stderr, "ERROR: clBuildProgram failure for %s\n",pName);
		fprintf(stderr, "Build Log for %s_program:\n%s\n", pName, build_c );
		sclPrintErrorFlags( err ); 
	}

}

cl_kernel _sclCreateKernel( sclSoft software ) {
	cl_kernel kernel;
	cl_int err;

	kernel = clCreateKernel( software.program, software.kernelName, &err );
	if ( err != CL_SUCCESS ) {
		printf( "Error on createKernel %s ", software.kernelName );
		fprintf(stderr, "Error on createKernel %s ", software.kernelName );
		sclPrintErrorFlags( err );
	}

	return kernel;
}



void sclEnqueueKernel( sclHard hardware, sclSoft software) {

	cl_int err;

	err = clEnqueueNDRangeKernel( hardware.queue, software.kernel, 3, NULL, software.global_size, software.local_size, 0, NULL, NULL );
	if ( err != CL_SUCCESS ) {
		printf( "\nError on EnqueueKernel %s", software.kernelName );
		fprintf(stderr, "\nError on EnqueueKernel %s", software.kernelName );
		sclPrintErrorFlags(err); 
	}
		
}


cl_event sclEnqueueKernelEvent( sclHard hardware, sclSoft software) {

	cl_event myEvent;
	cl_int err;

	err = clEnqueueNDRangeKernel( hardware.queue, software.kernel, 3, NULL, software.global_size, software.local_size, 0, NULL, &myEvent );
	if ( err != CL_SUCCESS ) {
		printf( "\nError on EnqueueKernel %s", software.kernelName );
		fprintf(stderr, "\nError on EnqueueKernel %s", software.kernelName );
		sclPrintErrorFlags(err); 
	}

	return myEvent;
		
}


double ProfilesclEnqueueKernel( sclHard hardware, sclSoft software) {
	cl_event myEvent;	
	cl_int err;
	cl_ulong time_start;
	cl_ulong time_end;
	double ms = 0.0;

	err = clEnqueueNDRangeKernel( hardware.queue, software.kernel, 3, NULL, software.global_size, software.local_size, 0, NULL, &myEvent );
	if ( err != CL_SUCCESS ) {
		printf( "\nError on EnqueueKernel %s", software.kernelName );
		fprintf(stderr, "\nError on EnqueueKernel %s", software.kernelName );
		sclPrintErrorFlags(err); 
	}

	clWaitForEvents(1, &myEvent);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	ms = (time_end-time_start) / 1000000.0;
	//printf("%0.3f ms %s kernel time\n", ms, software.kernelName);

	clReleaseEvent(myEvent);

	return ms;		
}


double ProfilesclEnqueueKernelNS( sclHard hardware, sclSoft software) {
	cl_event myEvent;	
	cl_int err;
	cl_ulong time_start;
	cl_ulong time_end;
	double ns = 0.0;

	err = clEnqueueNDRangeKernel( hardware.queue, software.kernel, 3, NULL, software.global_size, software.local_size, 0, NULL, &myEvent );
	if ( err != CL_SUCCESS ) {
		printf( "\nError on EnqueueKernel %s", software.kernelName );
		fprintf(stderr, "\nError on EnqueueKernel %s", software.kernelName );
		sclPrintErrorFlags(err); 
	}

	clWaitForEvents(1, &myEvent);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	ns = time_end-time_start;

	clReleaseEvent(myEvent);

	return ns;		
}



void sclSetGlobalSize( sclSoft & software, uint64_t size ) {

	uint64_t newsize = (size / software.local_size[0]) * software.local_size[0];

	if(newsize < size || newsize == 0){
		newsize += software.local_size[0];
	}

	software.global_size[0] = newsize;

}


void sclSetGlobalSizeExact( sclSoft & software, uint64_t size ) {

	software.global_size[0] = size;

}


void sclReleaseClSoft( sclSoft soft ) {
	clReleaseKernel( soft.kernel );
	clReleaseProgram( soft.program );
}


void sclReleaseClHard( sclHard hardware ){
	clReleaseCommandQueue( hardware.queue );
	clReleaseContext( hardware.context );
}


void sclReleaseMemObject( cl_mem object ) {
	cl_int err;

	if(object != NULL){
		err = clReleaseMemObject( object );
		if ( err != CL_SUCCESS ) {
			printf( "\nError on sclReleaseMemObject" );
			fprintf(stderr, "\nError on sclReleaseMemObject" );
			sclPrintErrorFlags(err); 
		}	
	}
}


int _sclGetMaxComputeUnits( cl_device_id device ) {
	
	cl_uint nCompU;

	clGetDeviceInfo( device, CL_DEVICE_MAX_COMPUTE_UNITS, 4, (void *)&nCompU, NULL );

	return (int)nCompU;	

}

unsigned long int _sclGetMaxMemAllocSize( cl_device_id device ){

	cl_ulong mem;

 	clGetDeviceInfo( device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, 8, (void *)&mem, NULL );

	return (unsigned long int)mem;	

}


unsigned long int _sclGetMaxGlobalMemSize( cl_device_id device ){

	cl_ulong mem;

 	clGetDeviceInfo( device, CL_DEVICE_GLOBAL_MEM_SIZE, 8, (void *)&mem, NULL );

	return (unsigned long int)mem;	

}


sclSoft sclGetCLSoftware( const char* source, const char* name, sclHard hardware, const char * options ){

	sclSoft software;

	sprintf( software.kernelName, "%s", name);
	
	/* Create program objects from source
	 ########################################################### */
	software.program = _sclCreateProgram( source, hardware.context );
	/* ########################################################### */
	
	/* Build the program (compile it)
   	 ############################################ */
   	_sclBuildProgram( software.program, hardware.device, name, options );
   	/* ############################################ */
   	
   	/* Create the kernel object
	 ########################################################################## */
	software.kernel = _sclCreateKernel( software );
	/* ########################################################################## */


	cl_int err;
	size_t workgroupsize;

	err = clGetKernelWorkGroupInfo( software.kernel, hardware.device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &workgroupsize, NULL);

	if ( err != CL_SUCCESS ) {
		printf( "\nError getting kernel workgroup size\n");
		fprintf(stderr, "Error getting kernel workgroup size\n");
		sclPrintErrorFlags(err); 
	}

	software.local_size[0] = workgroupsize;

	return software;
	
}

sclSoft sclGetCLSoftwareWithCommon( const char* common, const char* source, const char* name, sclHard hardware, const char * options ){

	sclSoft software;

	sprintf( software.kernelName, "%s", name);
	
	size_t total_len = strlen(common) + strlen(source) + 1;
	char *combined_src = (char *)malloc(total_len);
	strcpy(combined_src, common);
	strcat(combined_src, source);
	
	/* Create program objects from source
	 ########################################################### */
	software.program = _sclCreateProgram( (const char *)combined_src, hardware.context );
	/* ########################################################### */
	
        free(combined_src);	
	
	/* Build the program (compile it)
   	 ############################################ */
   	_sclBuildProgram( software.program, hardware.device, name, options );
   	/* ############################################ */
   	
   	/* Create the kernel object
	 ########################################################################## */
	software.kernel = _sclCreateKernel( software );
	/* ########################################################################## */


	cl_int err;
	size_t workgroupsize;

	err = clGetKernelWorkGroupInfo( software.kernel, hardware.device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &workgroupsize, NULL);

	if ( err != CL_SUCCESS ) {
		printf( "\nError getting kernel workgroup size\n");
		fprintf(stderr, "Error getting kernel workgroup size\n");
		sclPrintErrorFlags(err); 
	}

	software.local_size[0] = workgroupsize;

	return software;
	
}

void sclWrite( sclHard hardware, size_t size, cl_mem buffer, void* hostPointer ) {

	cl_int err;

	err = clEnqueueWriteBuffer( hardware.queue, buffer, CL_TRUE, 0, size, hostPointer, 0, NULL, NULL );
	if ( err != CL_SUCCESS ) { 
		printf( "\nclWrite Error\n" );
		fprintf(stderr, "\nclWrite Error\n" );
		sclPrintErrorFlags( err );
	}   

}

void sclWriteNB( sclHard hardware, size_t size, cl_mem buffer, void* hostPointer ) {

	cl_int err;

	err = clEnqueueWriteBuffer( hardware.queue, buffer, CL_FALSE, 0, size, hostPointer, 0, NULL, NULL );
	if ( err != CL_SUCCESS ) { 
		printf( "\nclWrite Error\n" );
		fprintf(stderr, "\nclWrite Error\n" );
		sclPrintErrorFlags( err );
	}   

}

void sclRead( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer ) {

	cl_int err;

	err = clEnqueueReadBuffer( hardware.queue, buffer, CL_TRUE, 0, size, hostPointer, 0, NULL, NULL );
	if ( err != CL_SUCCESS ) {
		printf( "\nclRead Error\n" );
		fprintf(stderr, "\nclRead Error\n" );
		sclPrintErrorFlags( err );
       	}

}

void sclReadNB( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer ) {

	cl_int err;

	err = clEnqueueReadBuffer( hardware.queue, buffer, CL_FALSE, 0, size, hostPointer, 0, NULL, NULL );
	if ( err != CL_SUCCESS ) {
		printf( "\nclRead Error\n" );
		fprintf(stderr, "\nclRead Error\n" );
		sclPrintErrorFlags( err );
       	}

}

cl_event sclReadNBEvent( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer ) {

	cl_int err;
	cl_event myEvent;	

	err = clEnqueueReadBuffer( hardware.queue, buffer, CL_FALSE, 0, size, hostPointer, 0, NULL, &myEvent );
	if ( err != CL_SUCCESS ) {
		printf( "\nclRead Error\n" );
		fprintf(stderr, "\nclRead Error\n" );
		sclPrintErrorFlags( err );
       	}

	return myEvent;
}

double ProfilesclRead( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer ) {
	cl_event myEvent;	
	cl_int err;
	cl_ulong time_start;
	cl_ulong time_end;
	double ns = 0.0;

	err = clEnqueueReadBuffer( hardware.queue, buffer, CL_TRUE, 0, size, hostPointer, 0, NULL, &myEvent );
	if ( err != CL_SUCCESS ) {
		printf( "\nclRead Error\n" );
		fprintf(stderr, "\nclRead Error\n" );
		sclPrintErrorFlags( err );
       	}

	clWaitForEvents(1, &myEvent);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	ns = time_end-time_start;

	clReleaseEvent(myEvent);

	return ns;		
}

cl_int sclFinish( sclHard hardware ){

	cl_int err;

	err = clFinish( hardware.queue );
	if ( err != CL_SUCCESS ) {
		printf( "\nError clFinish\n" );
		fprintf(stderr, "\nError clFinish\n" );
		sclPrintErrorFlags( err );
	}


	return err;

}


void sclSetKernelArg( sclSoft software, int argnum, size_t typeSize, void *argument ){

	cl_int err;

	err = clSetKernelArg( software.kernel, argnum, typeSize, argument );
	if ( err != CL_SUCCESS ) {	
		printf( "\nError clSetKernelArg number %d for %s\n", argnum, software.kernelName );
		fprintf(stderr, "\nError clSetKernelArg number %d for %s\n", argnum, software.kernelName );
		sclPrintErrorFlags( err );
	}


}


#ifdef __cplusplus
}
#endif
