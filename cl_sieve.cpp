/*
	PCWSieve
	Bryan Little, Jun 6 2025
	
	Search algorithm by
	Geoffrey Reynolds, 2009
	Ken Brazier, 2009
	https://github.com/Ken-g6/PSieve-CUDA/tree/redcl
	https://github.com/Ken-g6/PSieve-CUDA/tree/cw

	With contributions by
	Yves Gallot

*/

#include <unistd.h>
#include <math.h>

#include "boinc_api.h"
#include "boinc_opencl.h"
#include "simpleCL.h"

#include "clearn.h"
#include "clearresult.h"
#include "getsegprimes.h"
#include "sieve.h"
#include "sievecw.h"
#include "setup.h"
#include "check.h"

#include "factor_proth.h"
#include "verify_factor.h"
#include "putil.h"
#include "cl_sieve.h"

#define RESULTS_FILENAME "factors.txt"
#define STATE_FILENAME_A "PCWstateA.txt"
#define STATE_FILENAME_B "PCWstateB.txt"

void handle_trickle_up(searchData & sd)
{
	if(boinc_is_standalone()) return;

	uint64_t now = (uint64_t)time(NULL);

	if( (now-sd.last_trickle) > 86400 ){	// Once per day

		sd.last_trickle = now;

		double progress = boinc_get_fraction_done();
		double cpu;
		boinc_wu_cpu_time(cpu);
		APP_INIT_DATA init_data;
		boinc_get_init_data(init_data);
		double run = boinc_elapsed_time() + init_data.starting_elapsed_time;

		char msg[512];
		sprintf(msg, "<trickle_up>\n"
			    "   <progress>%lf</progress>\n"
			    "   <cputime>%lf</cputime>\n"
			    "   <runtime>%lf</runtime>\n"
			    "</trickle_up>\n",
			     progress, cpu, run  );
		char variety[64];
		sprintf(variety, "cwsieve_progress");
		boinc_send_trickle_up(variety, msg);
	}

}


FILE *my_fopen(const char * filename, const char * mode)
{
	char resolved_name[512];

	boinc_resolve_filename(filename,resolved_name,sizeof(resolved_name));
	return boinc_fopen(resolved_name,mode);

}


typedef struct {

	uint32_t range;
	uint32_t psize;
	uint32_t numgroups;

	cl_mem d_factorP = NULL;
	cl_mem d_factorKN = NULL;
	cl_mem d_factorcount = NULL;

	cl_mem d_flag = NULL;
	cl_mem d_checksum = NULL;

	cl_mem d_primes = NULL;
	cl_mem d_primecount = NULL;

	cl_mem d_Ps = NULL;
	cl_mem d_K = NULL;
	cl_mem d_lK = NULL;

	sclSoft sieve, clearn, clearresult, setup, check, getsegprimes;

}progData;


void cleanup( progData pd ){
	sclReleaseMemObject(pd.d_factorP);
	sclReleaseMemObject(pd.d_factorKN);
	sclReleaseMemObject(pd.d_factorcount);

	sclReleaseMemObject(pd.d_flag);
	sclReleaseMemObject(pd.d_checksum);

	sclReleaseMemObject(pd.d_primes);
	sclReleaseMemObject(pd.d_primecount);

	sclReleaseMemObject(pd.d_Ps);
	sclReleaseMemObject(pd.d_K);
	sclReleaseMemObject(pd.d_lK);

	sclReleaseClSoft(pd.clearn);
	sclReleaseClSoft(pd.clearresult);
        sclReleaseClSoft(pd.sieve);
        sclReleaseClSoft(pd.setup);
        sclReleaseClSoft(pd.check);
        sclReleaseClSoft(pd.getsegprimes);

}


void write_state( searchData & sd ){

	FILE *out;

        if (sd.write_state_a_next){
		if ((out = my_fopen(STATE_FILENAME_A,"w")) == NULL)
			fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_A);
	}
	else{
                if ((out = my_fopen(STATE_FILENAME_B,"w")) == NULL)
                        fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_B);
        }
	if (fprintf(out,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",sd.workunit,sd.p,sd.primecount,sd.checksum,sd.factorcount,sd.last_trickle) < 0){
		if (sd.write_state_a_next)
			fprintf(stderr,"Cannot write to %s !!! Continuing...\n",STATE_FILENAME_A);
		else
			fprintf(stderr,"Cannot write to %s !!! Continuing...\n",STATE_FILENAME_B);

		// Attempt to close, even though we failed to write
		fclose(out);
	}
	else{
		// If state file is closed OK, write to the other state file
		// next time around
		if (fclose(out) == 0) 
			sd.write_state_a_next = !sd.write_state_a_next; 
	}
}

/* Return 1 only if a valid checkpoint can be read.
   Attempts to read from both state files,
   uses the most recent one available.
 */
int read_state( searchData & sd ){

	FILE *in;
	bool good_state_a = true;
	bool good_state_b = true;
	uint64_t workunit_a, workunit_b;
	uint64_t current_a, current_b;
	uint64_t primecount_a, primecount_b;
	uint64_t checksum_a, checksum_b;
	uint64_t factorcount_a, factorcount_b;
	uint64_t trickle_a, trickle_b;

        // Attempt to read state file A
	if ((in = my_fopen(STATE_FILENAME_A,"r")) == NULL){
		good_state_a = false;
        }
	else if (fscanf(in,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",&workunit_a,&current_a,&primecount_a,&checksum_a,&factorcount_a,&trickle_a) != 6){
		fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_A);
		good_state_a = false;
	}
	else{
		fclose(in);
		/* Check that start stop match */
		if (workunit_a != sd.workunit){
			good_state_a = false;
		}
	}

        // Attempt to read state file B
        if ((in = my_fopen(STATE_FILENAME_B,"r")) == NULL){
                good_state_b = false;
        }
	else if (fscanf(in,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",&workunit_b,&current_b,&primecount_b,&checksum_b,&factorcount_b,&trickle_b) != 6){
                fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_B);
                good_state_b = false;
        }
        else{
                fclose(in);
		/* Check that start stop match */
		if (workunit_b != sd.workunit){
                        good_state_b = false;
                }
        }

        // If both state files are OK, check which is the most recent
	if (good_state_a && good_state_b)
	{
		if (current_a > current_b)
			good_state_b = false;
		else
			good_state_a = false;
	}

        // Use data from the most recent state file
	if (good_state_a && !good_state_b)
	{
		sd.p = current_a;
		sd.primecount = primecount_a;
		sd.checksum = checksum_a;
		sd.factorcount = factorcount_a;
		sd.last_trickle = trickle_a;
		sd.write_state_a_next = false;
		return 1;
	}
        if (good_state_b && !good_state_a)
        {
                sd.p = current_b;
		sd.primecount = primecount_b;
		sd.checksum = checksum_b;
		sd.factorcount = factorcount_b;
		sd.last_trickle = trickle_b;
		sd.write_state_a_next = true;
		return 1;
        }

	// If we got here, neither state file was good
	return 0;
}


void checkpoint( searchData & sd ){

	handle_trickle_up( sd );

	write_state( sd );

	if(boinc_is_standalone()){
		printf("Checkpoint, current p: %" PRIu64 "\n", sd.p);
	}

	boinc_checkpoint_completed();
}


// sleep CPU thread while waiting on the specified event to complete in the command queue
// using critical sections to prevent BOINC from shutting down the program while kernels are running on the GPU
void waitOnEvent(sclHard hardware, cl_event event){

	cl_int err;
	cl_int info;
#ifdef _WIN32
#else
	struct timespec sleep_time;
	sleep_time.tv_sec = 0;
	sleep_time.tv_nsec = 1000000;	// 1ms
#endif

	boinc_begin_critical_section();

	err = clFlush(hardware.queue);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clFlush\n" );
		fprintf(stderr, "ERROR: clFlush\n" );
		sclPrintErrorFlags( err );
       	}

	while(true){

#ifdef _WIN32
		Sleep(1);
#else
		nanosleep(&sleep_time,NULL);
#endif

		err = clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
		if ( err != CL_SUCCESS ) {
			printf( "ERROR: clGetEventInfo\n" );
			fprintf(stderr, "ERROR: clGetEventInfo\n" );
			sclPrintErrorFlags( err );
	       	}

		if(info == CL_COMPLETE){
			err = clReleaseEvent(event);
			if ( err != CL_SUCCESS ) {
				printf( "ERROR: clReleaseEvent\n" );
				fprintf(stderr, "ERROR: clReleaseEvent\n" );
				sclPrintErrorFlags( err );
		       	}

			boinc_end_critical_section();

			return;
		}
	}
}


// queue a marker and sleep CPU thread until marker has been reached in the command queue
void sleepCPU(sclHard hardware){

	cl_event kernelsDone;
	cl_int err;
	cl_int info;
#ifdef _WIN32
#else
	struct timespec sleep_time;
	sleep_time.tv_sec = 0;
	sleep_time.tv_nsec = 1000000;	// 1ms
#endif

	boinc_begin_critical_section();

	// OpenCL v2.0
/*
	err = clEnqueueMarkerWithWaitList( hardware.queue, 0, NULL, &kernelsDone);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clEnqueueMarkerWithWaitList\n");
		fprintf(stderr, "ERROR: clEnqueueMarkerWithWaitList\n");
		sclPrintErrorFlags(err); 
	}
*/
	err = clEnqueueMarker( hardware.queue, &kernelsDone);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clEnqueueMarker\n");
		fprintf(stderr, "ERROR: clEnqueueMarker\n");
		sclPrintErrorFlags(err); 
	}

	err = clFlush(hardware.queue);
	if ( err != CL_SUCCESS ) {
		printf( "ERROR: clFlush\n" );
		fprintf(stderr, "ERROR: clFlush\n" );
		sclPrintErrorFlags( err );
       	}

	while(true){

#ifdef _WIN32
		Sleep(1);
#else
		nanosleep(&sleep_time,NULL);
#endif

		err = clGetEventInfo(kernelsDone, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &info, NULL);
		if ( err != CL_SUCCESS ) {
			printf( "ERROR: clGetEventInfo\n" );
			fprintf(stderr, "ERROR: clGetEventInfo\n" );
			sclPrintErrorFlags( err );
	       	}

		if(info == CL_COMPLETE){
			err = clReleaseEvent(kernelsDone);
			if ( err != CL_SUCCESS ) {
				printf( "ERROR: clReleaseEvent\n" );
				fprintf(stderr, "ERROR: clReleaseEvent\n" );
				sclPrintErrorFlags( err );
		       	}

			boinc_end_critical_section();

			return;
		}
	}
}


// find mod 30 wheel index based on starting N
// this is used by gpu threads to iterate over the number line
void findWheelOffset(uint64_t & start, int32_t & index){

	int32_t wheel[8] = {4, 2, 4, 2, 4, 6, 2, 6};
	int32_t idx = -1;

	// find starting number using mod 6 wheel
	// N=(k*6)-1, N=(k*6)+1 ...
	// where k, k+1, k+2 ...
	uint64_t k = start / 6;
	int32_t i = 1;
	uint64_t N = (k * 6)-1;


	while( N < start || N % 5 == 0 ){
		if(i){
			i = 0;
			N += 2;
		}
		else{
			i = 1;
			N += 4;
		}
	}

	start = N;

	// find mod 30 wheel index by iterating with a mod 6 wheel until finding N divisible by 5
	// forward to find index
	while(idx < 0){

		if(i){
			N += 2;
			i = 0;
			if(N % 5 == 0){
				N -= 2;
				idx = 5;
			}

		}
		else{
			N += 4;
			i = 1;
			if(N % 5 == 0){
				N -= 4;
				idx = 7;
			}
		}
	}

	// reverse to find starting index
	while(N != start){
		--idx;
		if(idx < 0)idx=7;
		N -= wheel[idx];
	}


	index = idx;

}


void getResults( progData pd, searchData & sd, sclHard hardware ){

	uint64_t * h_checksum = (uint64_t *)malloc(pd.numgroups*sizeof(uint64_t));
	if( h_checksum == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy checksum and total prime count to host memory
	// blocking read
	sclRead(hardware, pd.numgroups*sizeof(uint64_t), pd.d_checksum, h_checksum);

	// index 0 is the gpu's total prime count
	sd.primecount += h_checksum[0];

	// sum block checksums
	for(uint32_t i=1; i<pd.numgroups; ++i){
		sd.checksum += h_checksum[i];
	}

	uint32_t * h_flag = (uint32_t *)malloc(sizeof(uint32_t));
	if( h_flag == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy checksum flag to host memory
	// blocking read
	sclRead(hardware, sizeof(uint32_t), pd.d_flag, h_flag);

	// flag set by gpu if there is an internal checksum error
	if(*h_flag > 0){
		fprintf(stderr,"error: gpu checksum failure\n");
		printf("error: gpu checksum failure\n");
		exit(EXIT_FAILURE);
	}

	uint32_t * h_primecount = (uint32_t *)malloc(2*sizeof(uint32_t));
	if( h_primecount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy prime count to host memory
	// blocking read
	sclRead(hardware, 2*sizeof(uint32_t), pd.d_primecount, h_primecount);

	// largest kernel prime count.  used to check array bounds
	if(h_primecount[1] > pd.psize){
		fprintf(stderr,"error: gpu prime array overflow\n");
		printf("error: gpu prime array overflow\n");
		exit(EXIT_FAILURE);
	}

	uint32_t * h_factorcount = (uint32_t *)malloc(sizeof(uint32_t));
	if( h_factorcount == NULL ){
		fprintf(stderr,"malloc error\n");
		exit(EXIT_FAILURE);
	}

	// copy factor cnt to host memory
	// blocking read
	sclRead(hardware, sizeof(uint32_t), pd.d_factorcount, h_factorcount);

//	printf("%u factors found on gpu.  verifying on cpu.\n",*h_factorcount);

	if(*h_factorcount > 0){

		if(*h_factorcount > sd.numresults){
			fprintf(stderr,"Error: number of results (%u) overflowed array.\n", *h_factorcount);
			exit(EXIT_FAILURE);
		}

		int64_t * h_factorP = (int64_t *)malloc(*h_factorcount * sizeof(int64_t));
		if( h_factorP == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}

		cl_uint2 * h_factorKN = (cl_uint2 *)malloc(*h_factorcount * sizeof(cl_uint2));
		if( h_factorKN == NULL ){
			fprintf(stderr,"malloc error\n");
			exit(EXIT_FAILURE);
		}

		// copy factors to host memory
		// blocking read
		sclRead(hardware, *h_factorcount * sizeof(int64_t), pd.d_factorP, h_factorP);
		sclRead(hardware, *h_factorcount * sizeof(cl_uint2), pd.d_factorKN, h_factorKN);

		// sort results by prime size if needed
		if(*h_factorcount > 1){
			for (uint32_t i = 0; i < *h_factorcount-1; i++){    
				for (uint32_t j = 0; j < *h_factorcount-i-1; j++){
					uint64_t a = (h_factorP[j]<0)?-h_factorP[j]:h_factorP[j];
					uint64_t b = (h_factorP[j+1]<0)?-h_factorP[j+1]:h_factorP[j+1];
					if (a > b){
						std::swap(h_factorP[j], h_factorP[j+1]);
						std::swap(h_factorKN[j], h_factorKN[j+1]);
					}
				}
			}
		}


		FILE * resfile = my_fopen(RESULTS_FILENAME,"a");

		if( resfile == NULL ){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}

		for(uint32_t m=0; m<*h_factorcount; ++m){

			int64_t sp;
			uint64_t p;
			uint32_t k;
			uint32_t n;
			int32_t c;

			// use the sign bit of P for the sign of the factor since its limited to 2^62
			sp = h_factorP[m];
			p = (sp < 0)?-sp:sp;
			k = h_factorKN[m].s0;
			n = h_factorKN[m].s1;
			c = (sp < 0)?-1:1;

			if(sd.cw){

				if(try_all_factors(k, n, c) == 0){	// check for a small prime factor of the number

					// check the factor actually divides the number
					if(verify_factor(p,k,n,c)){
						++sd.factorcount;
						if( fprintf( resfile, "%" PRIu64 " | %u*2^%u%+d\n",p,k,n,c) < 0 ){
							fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
							exit(EXIT_FAILURE);
						}						
						// add the factor to checksum
						sd.checksum += k;
						sd.checksum += n;
						(c == 1)?(++sd.checksum):(--sd.checksum);
					}
					else{
						printf("ERROR: GPU calculated invalid factor!\n");
						fprintf(stderr,"ERROR: GPU calculated invalid factor!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else{
				uint64_t b = k/sd.kstep;

				if(k == sd.kstep*b+sd.koffset) { // k is odd.

					if(try_all_factors(k, n, c) == 0 ){  // check for a small prime factor of the number

						// check the factor actually divides the number
						if(verify_factor(p,k,n,c)){
							++sd.factorcount;
							if( fprintf( resfile, "%" PRIu64 " | %u*2^%u%+d\n",p,k,n,c) < 0 ){
								fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
								exit(EXIT_FAILURE);
							}								
							// add the factor to checksum
							sd.checksum += k;
							sd.checksum += n;
							(c == 1)?(++sd.checksum):(--sd.checksum);
						}
						else{
							printf("ERROR: GPU calculated invalid factor!\n");
							fprintf(stderr,"ERROR: GPU calculated invalid factor!\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}

		}

		fclose(resfile);

		free(h_factorP);
		free(h_factorKN);
	}

	free(h_flag);
	free(h_factorcount);
	free(h_checksum);
	free(h_primecount);

}



// find the log base 2 of a number.
int lg2(uint64_t v) {

#ifdef __GNUC__
	return 63-__builtin_clzll(v);
#else
	int r = 0; // r will be lg(v)
	while (v >>= 1)r++;
	return r;
#endif

}


void setupSearch(searchData & sd){

	sd.p = sd.pmin;

	if(sd.pmin == 0 || sd.pmax == 0){
		printf("-p and -P arguments are required\n");
		fprintf(stderr, "-p and -P arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if(sd.nmin == 0 || sd.nmax == 0){
		printf("-n and -N arguments are required\n");
		fprintf(stderr, "-n and -N arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if (sd.nmin > sd.nmax){
		printf("nmin <= nmax is required\n");
		fprintf(stderr, "nmin <= nmax is required\n");
		exit(EXIT_FAILURE);
	}

	if(sd.cw){

		if(sd.nmax >= sd.pmin){
			printf("nmax < pmin is required\n");
			fprintf(stderr, "nmax < pmin is required\n");
			exit(EXIT_FAILURE);
		}

		sd.kmax = sd.nmax;
		sd.kmin = sd.nmin;
	}
	else{

		if(sd.kmax == 0){
			printf("-K argument is required\n");
			fprintf(stderr, "-K argument is required\n");
			exit(EXIT_FAILURE);
		}

		if(sd.kmin > sd.kmax){
			printf("kmin <= kmax is required\n");
			fprintf(stderr, "kmin <= kmax is required\n");
			exit(EXIT_FAILURE);
		}

		if(sd.kmax >= sd.pmin){
			printf("kmax < pmin is required\n");
			fprintf(stderr, "kmax < pmin is required\n");
			exit(EXIT_FAILURE);
		}

		uint32_t b0 = 0, b1 = 0;
		b0 = sd.kmin/sd.kstep;
		b1 = sd.kmax/sd.kstep;
		sd.kmin = b0*sd.kstep+sd.koffset;
		sd.kmax = b1*sd.kstep+sd.koffset;
	}


	for (sd.nstep = 1; ( (uint64_t)(sd.kmax) << sd.nstep ) < sd.pmin; sd.nstep++);

	if((((uint64_t)1) << (64-sd.nstep)) > sd.pmin) {

		uint64_t pmin_1 = (((uint64_t)1) << (64-sd.nstep));

		printf("Error: pmin is not large enough (or nmax is close to nmin).\n");
		fprintf(stderr, "Error: pmin is not large enough (or nmax is close to nmin).\n");

		sd.pmin = sd.kmax + 1;
		for (sd.nstep = 1; ( (uint64_t)(sd.kmax) << sd.nstep ) < sd.pmin; sd.nstep++);

		while((((uint64_t)1) << (64-sd.nstep)) > sd.pmin) {
			sd.pmin *= 2;
			sd.nstep++;
		}
		if(pmin_1 < sd.pmin){
			sd.pmin = pmin_1;
		}

		printf("This program will work by the time pmin == %" PRIu64 ".\n", sd.pmin);
		fprintf(stderr, "This program will work by the time pmin == %" PRIu64 ".\n", sd.pmin);

		exit(EXIT_FAILURE);
	}

	if (sd.nstep > (sd.nmax-sd.nmin+1))
		sd.nstep = (sd.nmax-sd.nmin+1);

	// For TPS, decrease the ld_nstep by one to allow overlap, checking both + and -
	sd.nstep--;

	// Use the 32-step algorithm where useful.
	if(sd.nstep >= 32 && (((uint64_t)1) << 32) <= sd.pmin) {
		sd.nstep = 32;
	}

	// N's to search each time a kernel is run
	if(sd.compute){
		sd.kernel_nstep = sd.nstep * 15000;
	}
	else{
		sd.kernel_nstep = sd.nstep * 3000;
	}

	// search twin, decrease by one to allow overlap, checking both + and -
	sd.nmin--;

	sd.bbits = lg2(sd.nmin);

	if(sd.bbits < 6) {
		printf("Error: nmin too small at %d (must be at least 65).\n", sd.nmin+1);
		fprintf(stderr, "Error: nmin too small at %d (must be at least 65).\n", sd.nmin+1);
		exit(EXIT_FAILURE);
	}

	// r = 2^-i * 2^64 (mod N), something that can be done in a uint64_t!
	// If i is large (and it should be at least >= 32), there's a very good chance no mod is needed!
	sd.r0 = ((uint64_t)1) << (64-(sd.nmin >> (sd.bbits-5)));

	sd.bbits = sd.bbits-6;

	sd.mont_nstep = 64-sd.nstep;

	// data for checksum
	uint32_t maxn;

	maxn = ( (sd.nmax-sd.nmin) / sd.nstep) * sd.nstep;
	maxn += sd.nmin;

	if( maxn < sd.nmax ){
		maxn += sd.nstep;
	}

	int bbits1 = lg2(maxn) - 5;
	sd.r1 = ((uint64_t)1) << (64-(maxn >> bbits1));
	--bbits1;
	sd.bbits1 = bbits1;
	sd.lastN = maxn;

	// for checkpoints
	sd.workunit = sd.pmin + sd.pmax + (uint64_t)sd.nmin + (uint64_t)sd.nmax + (uint64_t)sd.kmin + (uint64_t)sd.kmax;
	
	// increase result buffer at low P range
	// it's still possible to overflow this with a fast GPU and large search range
	if(sd.pmin < 0xFFFFFFFF){
		sd.numresults = 30000000;
	}
	else{
		sd.numresults = 1000000;
	}


}



void profileGPU(progData & pd, searchData sd, sclHard hardware, int debuginfo ){

	// calculate approximate chunk size based on gpu's compute units
	cl_int err = 0;

	uint64_t calc_range = sd.computeunits * 750000;

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	uint64_t estimated = calc_range;

	uint64_t prof_start = sd.p;

	uint64_t prof_stop = prof_start + calc_range;

	sclSetGlobalSize( pd.getsegprimes, (calc_range/60)+1 );

	// get a count of primes in the gpu worksize
	uint64_t prof_range_primes = (prof_stop / log(prof_stop)) - (prof_start / log(prof_start));

	// calculate prime array size based on result
	uint64_t prof_mem_size = (uint64_t)(1.5 * (double)prof_range_primes);

	// kernels use uint for global id
	if(prof_mem_size > UINT32_MAX){
		fprintf(stderr, "ERROR: prof_mem_size too large.\n");
                printf( "ERROR: prof_mem_size too large.\n" );
		exit(EXIT_FAILURE);
	}

	// allocate temporary gpu prime array for profiling
	cl_mem d_profileprime = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, prof_mem_size*sizeof(uint64_t), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
	        printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

	int32_t prof_wheelidx;
	uint64_t prof_kernel_start = prof_start;

	findWheelOffset(prof_kernel_start, prof_wheelidx);

	// set static args
	sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &prof_kernel_start);
	sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &prof_stop);
	sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &prof_wheelidx);
	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &d_profileprime);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &pd.d_primecount);

	// zero prime count
	sclEnqueueKernel(hardware, pd.clearn);

	// Benchmark the GPU
	double kernel_ms = ProfilesclEnqueueKernel(hardware, pd.getsegprimes);

	// target is 10ms for prime generator kernel
	double prof_multi = 10.0 / kernel_ms;

	// update chunk size based on the profile
	calc_range = (uint64_t)( (double)calc_range * prof_multi );

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	if(debuginfo){
		printf("Kernel profile: %0.3f ms. Estimated / Actual worksize: %" PRIu64 " / %" PRIu64 "\n",kernel_ms,estimated,calc_range);
	}

	// get a count of primes in the new gpu worksize
	prof_stop = prof_start+calc_range;
	uint64_t range_primes = (prof_stop / log(prof_stop)) - (prof_start / log(prof_start));

	// calculate prime array size based on result
	uint64_t mem_size = (uint64_t)( 1.5 * (double)range_primes );

	if(mem_size > UINT32_MAX){
		fprintf(stderr, "ERROR: mem_size too large.\n");
                printf( "ERROR: mem_size too large.\n" );
		exit(EXIT_FAILURE);
	}

	pd.range = calc_range;
	pd.psize = mem_size;

	// free temporary array
	sclReleaseMemObject(d_profileprime);

}



void cl_sieve( sclHard hardware, searchData & sd ){

	progData pd;
	bool profile = true;
	bool debuginfo = false;
	time_t boinc_last, boinc_curr;
	time_t ckpt_curr, ckpt_last;
	cl_int err = 0;

	sieve_small_primes(11);

	// setup kernel parameters
	setupSearch(sd);

	fprintf(stderr, "Starting sieve at p: %" PRIu64 " n: %u k: %u\nStopping sieve at P: %" PRIu64 " N: %u K: %u\n", sd.pmin, sd.nmin+1, sd.kmin, sd.pmax, sd.nmax, sd.kmax);
	if(boinc_is_standalone()){
		printf("Starting sieve at p: %" PRIu64 " n: %u k: %u\nStopping sieve at P: %" PRIu64 " N: %u K: %u\n", sd.pmin, sd.nmin+1, sd.kmin, sd.pmax, sd.nmax, sd.kmax);
	}


	// device arrays
	pd.d_primecount = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, 2*sizeof(cl_uint), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_flag = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sizeof(cl_uint), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_factorP = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numresults*sizeof(cl_long), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_factorKN = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numresults*sizeof(cl_uint2), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_factorcount = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sizeof(cl_uint), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}


	if(sd.cw){
		if(sd.nstep == 32){
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecw32",hardware, 1, debuginfo);
		}
		else if(sd.nstep < 32){
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecwsm",hardware, 1, debuginfo);
		}
		else{
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecw",hardware, 1, debuginfo);
		}
	}
	else{
		if(sd.nstep == 32){
			pd.sieve = sclGetCLSoftware(sieve_cl,"sieve32",hardware, 1, debuginfo);
		}
		else if(sd.nstep < 32){
			pd.sieve = sclGetCLSoftware(sieve_cl,"sievesm",hardware, 1, debuginfo);
		}
		else{
			pd.sieve = sclGetCLSoftware(sieve_cl,"sieve",hardware, 1, debuginfo);
		}
	}

        pd.clearn = sclGetCLSoftware(clearn_cl,"clearn",hardware, 1, debuginfo);

        pd.clearresult = sclGetCLSoftware(clearresult_cl,"clearresult",hardware, 1, debuginfo);

        pd.setup = sclGetCLSoftware(setup_cl,"setup",hardware, 1, debuginfo);

        pd.check = sclGetCLSoftware(check_cl,"check",hardware, 1, debuginfo);

        pd.getsegprimes = sclGetCLSoftware(getsegprimes_cl,"getsegprimes",hardware, 1, debuginfo);


	// kernels have __attribute__ ((reqd_work_group_size(256, 1, 1)))
	// it's still possible the CL complier picked a different size
	if(pd.getsegprimes.local_size[0] != 256){
		pd.getsegprimes.local_size[0] = 256;
		fprintf(stderr, "Set getsegprimes kernel local size to 256\n");
	}
	if(pd.check.local_size[0] != 256){
		pd.check.local_size[0] = 256;
		fprintf(stderr, "Set check kernel local size to 256\n");
	}


	if( sd.test ){
		// clear result file
		FILE * temp_file = my_fopen(RESULTS_FILENAME,"w");
		if (temp_file == NULL){
			fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
		fclose(temp_file);
	}
	else{
		// Resume from checkpoint if there is one
		if( read_state( sd ) ){
			if(boinc_is_standalone()){
				printf("Resuming search from checkpoint. Current p: %" PRIu64 "\n", sd.p);
			}
			fprintf(stderr,"Resuming search from checkpoint. Current p: %" PRIu64 "\n", sd.p);

			//trying to resume a finished workunit
			if( sd.p == sd.pmax ){
				if(boinc_is_standalone()){
					printf("Workunit complete.\n");
				}
				fprintf(stderr,"Workunit complete.\n");
				boinc_finish(EXIT_SUCCESS);
			}
		}
		// starting from beginning
		else{
			// clear result file
			FILE * temp_file = my_fopen(RESULTS_FILENAME,"w");
			if (temp_file == NULL){
				fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
				exit(EXIT_FAILURE);
			}
			fclose(temp_file);

			// setup boinc trickle up
			sd.last_trickle = (uint64_t)time(NULL);
		}
	}

	// kernel used in profileGPU, setup arg
	sclSetKernelArg(pd.clearn, 0, sizeof(cl_mem), &pd.d_primecount);
	sclSetGlobalSize( pd.clearn, 64 );

	profileGPU(pd,sd,hardware,debuginfo);

	// number of gpu workgroups, used to size the checksum array on gpu
	pd.numgroups = (pd.psize / pd.check.local_size[0]) + 2;

	sclSetGlobalSize( pd.getsegprimes, (pd.range/60)+1 );
	sclSetGlobalSize( pd.setup, pd.psize );
	sclSetGlobalSize( pd.sieve, pd.psize );
	sclSetGlobalSize( pd.check, pd.psize );
	sclSetGlobalSize( pd.clearresult, pd.numgroups );

	// allocate gpu P, Ps, K, lastK arrays
	pd.d_primes = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, pd.psize*sizeof(cl_ulong), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_Ps = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, pd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_K = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, pd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_lK = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, pd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_checksum = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, pd.numgroups*sizeof(cl_ulong), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}


	// set static kernel args
	sclSetKernelArg(pd.clearresult, 0, sizeof(cl_mem), &pd.d_flag);
	sclSetKernelArg(pd.clearresult, 1, sizeof(cl_mem), &pd.d_factorcount);
	sclSetKernelArg(pd.clearresult, 2, sizeof(cl_mem), &pd.d_checksum);
	sclSetKernelArg(pd.clearresult, 3, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.clearresult, 4, sizeof(uint32_t), &pd.numgroups);

	////////////////////////
	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &pd.d_primecount);
	////////////////////////

	sclSetKernelArg(pd.setup, 0, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.setup, 1, sizeof(cl_mem), &pd.d_Ps);
	sclSetKernelArg(pd.setup, 2, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.setup, 3, sizeof(cl_mem), &pd.d_lK);
	sclSetKernelArg(pd.setup, 4, sizeof(uint64_t), &sd.r0);
	sclSetKernelArg(pd.setup, 5, sizeof(int32_t), &sd.bbits);
	sclSetKernelArg(pd.setup, 6, sizeof(uint32_t), &sd.nmin);
	sclSetKernelArg(pd.setup, 7, sizeof(uint64_t), &sd.r1);
	sclSetKernelArg(pd.setup, 8, sizeof(int32_t), &sd.bbits1);
	sclSetKernelArg(pd.setup, 9, sizeof(uint32_t), &sd.lastN);
	sclSetKernelArg(pd.setup, 10, sizeof(cl_mem), &pd.d_primecount);
	////////////////////////

	sclSetKernelArg(pd.sieve, 0, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.sieve, 1, sizeof(cl_mem), &pd.d_Ps);
	sclSetKernelArg(pd.sieve, 2, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.sieve, 3, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.sieve, 4, sizeof(cl_mem), &pd.d_factorKN);
	sclSetKernelArg(pd.sieve, 5, sizeof(cl_mem), &pd.d_factorP);
	sclSetKernelArg(pd.sieve, 6, sizeof(cl_mem), &pd.d_factorcount);

	sclSetKernelArg(pd.sieve, 8, sizeof(uint32_t), &sd.nstep);
	sclSetKernelArg(pd.sieve, 9, sizeof(uint32_t), &sd.kernel_nstep);
	sclSetKernelArg(pd.sieve, 10, sizeof(uint32_t), &sd.mont_nstep);
	sclSetKernelArg(pd.sieve, 11, sizeof(uint32_t), &sd.nmax);
	sclSetKernelArg(pd.sieve, 12, sizeof(uint32_t), &sd.kmin);
	sclSetKernelArg(pd.sieve, 13, sizeof(uint32_t), &sd.kmax);
	////////////////////////

	sclSetKernelArg(pd.check, 0, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.check, 1, sizeof(cl_mem), &pd.d_lK);
	sclSetKernelArg(pd.check, 2, sizeof(cl_mem), &pd.d_flag);
	sclSetKernelArg(pd.check, 3, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.check, 4, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.check, 5, sizeof(cl_mem), &pd.d_checksum);
	sclSetKernelArg(pd.check, 6, sizeof(uint32_t), &pd.numgroups);
	////////////////////////


	fprintf(stderr,"Starting search...\n");
	if(boinc_is_standalone()){
		printf("Starting search...\n");
	}

	time(&boinc_last);
	time(&ckpt_last);

	printf("nstep: %u\n",sd.nstep);

	// clear results, checksum, total prime counts
	sclEnqueueKernel(hardware, pd.clearresult);

	time_t totals, totalf;
	if(boinc_is_standalone()){
		time(&totals);
	}

	// main search loop
	for(uint64_t stop; sd.p < sd.pmax; sd.p = stop){

		// clear prime count
		sclEnqueueKernel(hardware, pd.clearn);

		stop = sd.p + pd.range;
		if(stop > sd.pmax) stop = sd.pmax;

		// update BOINC fraction done every 2 sec
		time(&boinc_curr);
		if( ((int)boinc_curr - (int)boinc_last) > 1 ){
    			double fd = (double)(sd.p-sd.pmin)/(double)(sd.pmax-sd.pmin);
			boinc_fraction_done(fd);
			if(boinc_is_standalone()) printf("Tests done: %.1f%%\n",fd*100.0);
			boinc_last = boinc_curr;
		}

		// 1 minute checkpoint
		time(&ckpt_curr);
		if( ((int)ckpt_curr - (int)ckpt_last) > 60 ){
			sleepCPU(hardware);
			boinc_begin_critical_section();
			getResults(pd, sd, hardware);
			checkpoint(sd);
			boinc_end_critical_section();
			ckpt_last = ckpt_curr;
			// clear result arrays
			sclEnqueueKernel(hardware, pd.clearresult);
		}

		// get primes
		int32_t wheelidx;
		uint64_t kernel_start = sd.p;
		findWheelOffset(kernel_start, wheelidx);

		sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &kernel_start);
		sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &stop);
		sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &wheelidx);
		cl_event launchEvent = sclEnqueueKernelEvent(hardware, pd.getsegprimes);

		// setup Ps, K kernel
		sclEnqueueKernel(hardware, pd.setup);

		uint32_t nstart = sd.nmin;

		// profile gpu sieve kernel time once, at program start.  adjust work size to target kernel runtime.
		if(profile){
			sclSetKernelArg(pd.sieve, 7, sizeof(uint32_t), &nstart);
			double kernel_ms = ProfilesclEnqueueKernel(hardware, pd.sieve);
			nstart += sd.kernel_nstep;
			double multi = (sd.compute)?(50.0 / kernel_ms):(10.0 / kernel_ms);	// target kernel time 50ms or 10ms
			uint32_t new_knstep = (uint32_t)((double)sd.kernel_nstep * multi);
			// make sure it's a multiple of nstep
			new_knstep = (new_knstep / sd.nstep) * sd.nstep;
			if(debuginfo) printf("old kns %u, new kns %u\n",sd.kernel_nstep,new_knstep);
			sd.kernel_nstep = new_knstep;
			sclSetKernelArg(pd.sieve, 9, sizeof(uint32_t), &sd.kernel_nstep);
			profile = false;
		}

		// sieve kernel, loop to nmax
		for(; nstart <= sd.nmax; nstart += sd.kernel_nstep){
			sclSetKernelArg(pd.sieve, 7, sizeof(uint32_t), &nstart);
			sclEnqueueKernel(hardware, pd.sieve);
//			float kernel_ms = ProfilesclEnqueueKernel(hardware, pd.sieve);
//			printf("sieve kernel time %0.2fms\n",kernel_ms);
		}

		// validate checksum kernel
		sclEnqueueKernel(hardware, pd.check);

		// limit cl queue depth and sleep cpu
		waitOnEvent(hardware, launchEvent);

	}


	// final checkpoint
	sleepCPU(hardware);
	boinc_begin_critical_section();
	sd.p = sd.pmax;
	boinc_fraction_done(1.0);
	if(boinc_is_standalone()) printf("Tests done: %.1f%%\n",100.0);
	getResults(pd, sd, hardware);
	checkpoint(sd);

	
	// print checksum
	FILE * resfile = my_fopen(RESULTS_FILENAME,"a");

	if(resfile == NULL){
		fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
		exit(EXIT_FAILURE);
	}

	if(sd.factorcount){
		if( fprintf( resfile, "%016" PRIX64 "\n", sd.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}
	else{
		if( fprintf( resfile, "no factors\n%016" PRIX64 "\n", sd.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}

	fclose(resfile);
	
	boinc_end_critical_section();


	fprintf(stderr,"Search complete.\nfactors %" PRIu64 ", prime count %" PRIu64 "\n", sd.factorcount, sd.primecount);

	if(boinc_is_standalone()){
		time(&totalf);
		printf("Search finished in %d sec.\n", (int)totalf - (int)totals);
		printf("factors %" PRIu64 ", prime count %" PRIu64 ", checksum %016" PRIX64 "\n", sd.factorcount, sd.primecount, sd.checksum);
	}


	cleanup(pd);

	small_primes_free();
}


void run_test( sclHard hardware, searchData & sd ){

	int goodtest = 0;

	printf("Beginning self test of 4 ranges.\n");

//	-p 25636026e6 -P 25636030e6 -n 10000000 -N 25000000 -c		nstep 19
	sd.pmin = 25636026000000;
	sd.pmax = 25636030000000;
	sd.nmin = 10000000;
	sd.nmax = 25000000;
	sd.kmin = 0;
	sd.kmax = 0;
	sd.cw = true;
	cl_sieve( hardware, sd );
	if( sd.factorcount == 2 && sd.primecount == 129869 && sd.checksum == 0x4544591DC69ACD83 ){
		printf("CW test case 1 passed.\n\n");
		fprintf(stderr,"CW test case 1 passed.\n");
		++goodtest;
	}
	else{
		printf("CW test case 1 failed.\n\n");
		fprintf(stderr,"CW test case 1 failed.\n");
	}
	sd.checksum = 0;
	sd.primecount = 0;
	sd.factorcount = 0;

//	-p 556439300e6 -P 556439440e6 -n 100 -N 100000 -c		nstep 32
	sd.pmin = 556439300000000;
	sd.pmax = 556439440000000;
	sd.nmin = 100;
	sd.nmax = 100000;
	sd.kmin = 0;
	sd.kmax = 0;
	sd.cw = true;
	cl_sieve( hardware, sd );
	if( sd.factorcount == 1 && sd.primecount == 4123452 && sd.checksum == 0x8FEC30979896A3C0 ){
		printf("CW test case 2 passed.\n\n");
		fprintf(stderr,"CW test case 2 passed.\n");
		++goodtest;
	}
	else{
		printf("CW test case 2 failed.\n\n");
		fprintf(stderr,"CW test case 2 failed.\n");
	}
	sd.checksum = 0;
	sd.primecount = 0;
	sd.factorcount = 0;


//	-p838338347800e6 -P838338347820e6 -k5 -K9999 -n6000000 -N9000000	nstep 32
	sd.pmin = 838338347800000000;
	sd.pmax = 838338347820000000;
	sd.nmin = 6000000;
	sd.nmax = 9000000;
	sd.kmin = 5;
	sd.kmax = 9999;
	sd.cw = false;
	cl_sieve( hardware, sd );
	if( sd.factorcount == 1 && sd.primecount == 484024 && sd.checksum == 0xA7DC855BCB311759 ){
		printf("test case 3 passed.\n\n");
		fprintf(stderr,"test case 3 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 3 failed.\n\n");
		fprintf(stderr,"test case 3 failed.\n");
	}
	sd.checksum = 0;
	sd.primecount = 0;
	sd.factorcount = 0;

//	-p42070000e6 -P42070050e6 -k 1201 -K 9999 -n 100 -N 2000000		nstep 31
	sd.pmin = 42070000000000;
	sd.pmax = 42070050000000;
	sd.nmin = 100;
	sd.nmax = 2000000;
	sd.kmin = 1201;
	sd.kmax = 9999;
	sd.cw = false;
	cl_sieve( hardware, sd );
	if( sd.factorcount == 70 && sd.primecount == 1592285 && sd.checksum == 0x727796B2D3677937 ){
		printf("test case 4 passed.\n\n");
		fprintf(stderr,"test case 4 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 4 failed.\n\n");
		fprintf(stderr,"test case 4 failed.\n");
	}



	if(goodtest == 4){
		printf("All test cases completed successfully!\n");
		fprintf(stderr, "All test cases completed successfully!\n");
	}
	else{
		printf("Self test FAILED!\n");
		fprintf(stderr, "Self test FAILED!\n");
	}

}


