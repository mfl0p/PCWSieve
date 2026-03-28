/* #######################################################################
    Copyright 2011 Oscar Amoros Huguet, Cristian Garcia Marin

    This file is part of SimpleOpenCL.

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

#define CL_TARGET_OPENCL_VERSION 110

#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include "cl.h"
#endif

#define WORKGROUP_X 64
#define WORKGROUP_Y 2

#ifndef _OCLUTILS_STRUCTS
typedef struct {
	cl_platform_id platform;
	cl_context context;
	cl_device_id device;
	cl_command_queue queue;
}sclHard;
typedef sclHard* ptsclHard;
typedef struct {
	cl_program program;
	cl_kernel kernel;
	char kernelName[98];

        size_t global_size[3] = {1,1,1};
        size_t local_size[3] = {1,1,1};

}sclSoft;
#define _OCLUTILS_STRUCTS
#endif

/* USER FUNCTIONS */

void sclGetBinary( sclSoft software );
void sclSetGlobalSize( sclSoft & software, uint64_t size );
void sclSetGlobalSizeExact( sclSoft & software, uint64_t size );

/* ####### Device memory allocation read and write  ####### */

//cl_mem 			sclMalloc( sclHard hardware, cl_int mode, size_t size );
void 			sclWrite( sclHard hardware, size_t size, cl_mem buffer, void* hostPointer );
void 			sclWriteNB( sclHard hardware, size_t size, cl_mem buffer, void* hostPointer );
void			sclReadNB( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer );
cl_event 		sclReadNBEvent( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer );
void			sclRead( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer );
double 			ProfilesclRead( sclHard hardware, size_t size, cl_mem buffer, void *hostPointer );

/* ######################################################## */

/* ####### inicialization of sclSoft structs  ############## */
sclSoft 		sclGetCLSoftware( const char* source, const char* name, sclHard hardware, const char * options );
sclSoft 		sclGetCLSoftwareWithCommon( const char* common, const char* source, const char* name, sclHard hardware, const char * options );
/* ######################################################## */

/* ####### Release and retain OpenCL objects ############## */

void 			sclReleaseClSoft( sclSoft soft );
void 			sclReleaseClHard( sclHard hard );
void			sclReleaseMemObject( cl_mem object );

/* ######################################################## */

/* ####### Debug functions ################################ */

void 			sclPrintErrorFlags( cl_int flag );

/* ######################################################## */

/* ####### Device execution ############################### */

void			sclEnqueueKernel( sclHard hardware, sclSoft software );
cl_event		sclEnqueueKernelEvent( sclHard hardware, sclSoft software );
double			ProfilesclEnqueueKernel( sclHard hardware, sclSoft software );
double			ProfilesclEnqueueKernelNS( sclHard hardware, sclSoft software );

/* ######################################################## */

/* ####### Queue management ############################### */

cl_int			sclFinish( sclHard hardware );

/* ######################################################## */

/* ####### Kernel argument setting ######################## */

void 			sclSetKernelArg( sclSoft software, int argnum, size_t typeSize, void *argument );

/* ######################################################## */

/* ####### Hardware init and selection #################### */

/* ######################################################## */

/* INTERNAL FUNCITONS */

/* ####### cl software management ######################### */
void 			_sclBuildProgram( cl_program program, cl_device_id devices, const char* pName, const char * options );
cl_kernel 		_sclCreateKernel( sclSoft software );
cl_program 		_sclCreateProgram( const char* program_source, cl_context context );
char* 			_sclLoadProgramSource( const char *filename );

/* ######################################################## */

/* ####### hardware management ############################ */

int 			_sclGetMaxComputeUnits( cl_device_id device );
unsigned long int 	_sclGetMaxMemAllocSize( cl_device_id device );
unsigned long int 	_sclGetMaxGlobalMemSize( cl_device_id device );


/* ######################################################## */

#ifdef __cplusplus
}
#endif
