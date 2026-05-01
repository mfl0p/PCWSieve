/*
	PCWSieve
	Bryan Little, Apr 2026
	
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
#include <ctime>
#include <cinttypes>

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
#define STATE_FILENAME_A "stateA.ckp"
#define STATE_FILENAME_B "stateB.ckp"

void handle_trickle_up(workStatus & st)
{
	if(boinc_is_standalone()) return;

	uint64_t now = (uint64_t)time(NULL);

	if( (now-st.last_trickle) > 86400 ){	// Once per day

		st.last_trickle = now;

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


FILE *my_fopen(const char * filename, const char * mode){
	char resolved_name[512];
	boinc_resolve_filename(filename,resolved_name,sizeof(resolved_name));
	return boinc_fopen(resolved_name,mode);
}

void cleanup( progData & pd ){
	sclReleaseMemObject(pd.d_factor);
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


void format_rate(double val, char *buf)
{
    const char *units[] = {"", "k", "M", "G", "T", "P"};
    int u = 0;

    while(val >= 1000.0 && u < 5){
        val /= 1000.0;
        u++;
    }

    sprintf(buf, "%6.2f %sp/s", val, units[u]);
}

void format_eta(double seconds, char *buf)
{
    int h = seconds / 3600;
    int m = ((int)seconds % 3600) / 60;
    int s = (int)seconds % 60;

    sprintf(buf, "%02d:%02d:%02d", h, m, s);
}

void print_progress(workStatus &st,
                    double *smooth_rate,
                    time_t start_time,
                    uint64_t run_start_p){
    const int bar_width = 40;

    uint64_t total_range = st.pmax - st.pmin;
    uint64_t done_total  = st.p - st.pmin;

    double progress = 0.0;
    if(total_range > 0)
        progress = (double)done_total / (double)total_range;

    if(progress < 0) progress = 0;
    if(progress > 1) progress = 1;

    int pos = (int)(progress * bar_width);

    time_t now = time(NULL);
    double elapsed = difftime(now, start_time);

    uint64_t done_this_run = 0;
    if(st.p > run_start_p)
        done_this_run = st.p - run_start_p;

    double inst_rate = (double)done_this_run / (elapsed > 0 ? elapsed : 1.0);

    /*
        Exponential smoothing.
    */
    if(*smooth_rate == 0.0)
        *smooth_rate = inst_rate;
    else
        *smooth_rate = 0.9 * (*smooth_rate) + 0.1 * inst_rate;

    uint64_t remaining = 0;
    if(st.pmax > st.p)
        remaining = st.pmax - st.p;

    double remain = (double)remaining / (*smooth_rate > 0.0 ? *smooth_rate : 1.0);

    char rate_str[32];
    char eta_str[32];

    format_rate(*smooth_rate, rate_str);
    format_eta(remain, eta_str);

    printf("\r[");

    for(int i = 0; i < bar_width; i++)
        putchar(i < pos ? '#' : '-');

    printf("] %6.2f%% | %s | ETA %s",
           progress * 100.0,
           rate_str,
           eta_str);

    fflush(stdout);
}

// using fast binary checkpoint files with checksum calculation
void write_state( workStatus & st, searchData & sd ){

	FILE * out;

	st.state_sum = st.pmin+st.pmax+st.p+st.checksum+st.primecount+st.factorcount+st.last_trickle+st.nmin+st.nmax+st.kmin+st.kmax+st.cw;

        if (sd.write_state_a_next){
		if ((out = my_fopen(STATE_FILENAME_A,"wb")) == NULL)
			fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_A);
	}
	else{
                if ((out = my_fopen(STATE_FILENAME_B,"wb")) == NULL)
                        fprintf(stderr,"Cannot open %s !!!\n",STATE_FILENAME_B);
        }

	if(out != NULL){

		if( fwrite(&st, sizeof(workStatus), 1, out) != 1 ){
			fprintf(stderr,"Cannot write checkpoint to file. Continuing...\n");
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
}

int read_state( workStatus & st, searchData & sd ){

	FILE * in;
	bool good_state_a = true;
	bool good_state_b = true;
	workStatus stat_a, stat_b;

        // Attempt to read state file A
	if ((in = my_fopen(STATE_FILENAME_A,"rb")) == NULL){
		good_state_a = false;
        }
	else{
		if( fread(&stat_a, sizeof(workStatus), 1, in) != 1 ){
			fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_A);
			printf("Cannot parse %s !!!\n",STATE_FILENAME_A);
			good_state_a = false;
		}
		else if(stat_a.pmin != st.pmin || stat_a.pmax != st.pmax || stat_a.nmin != st.nmin || stat_a.nmax != st.nmax ||
			stat_a.kmin != st.kmin || stat_a.kmax != st.kmax || stat_a.cw != st.cw){
			fprintf(stderr,"Invalid checkpoint file %s !!!\n",STATE_FILENAME_A);
			printf("Invalid checkpoint file %s !!!\n",STATE_FILENAME_A);
			good_state_a = false;
		}
		else{
			uint64_t state_sum = stat_a.pmin+stat_a.pmax+stat_a.p+stat_a.checksum+stat_a.primecount+stat_a.factorcount+
						stat_a.last_trickle+stat_a.nmin+stat_a.nmax+stat_a.kmin+stat_a.kmax+stat_a.cw;

			if(state_sum != stat_a.state_sum){
				fprintf(stderr,"Checksum error in %s !!!\n",STATE_FILENAME_A);
				printf("Checksum error in %s !!!\n",STATE_FILENAME_A);
				good_state_a = false;
			}
		}
		fclose(in);
	}

        // Attempt to read state file B
	if ((in = my_fopen(STATE_FILENAME_B,"rb")) == NULL){
		good_state_b = false;
        }
	else{
		if( fread(&stat_b, sizeof(workStatus), 1, in) != 1 ){
			fprintf(stderr,"Cannot parse %s !!!\n",STATE_FILENAME_B);
			printf("Cannot parse %s !!!\n",STATE_FILENAME_B);
			good_state_b = false;
		}
		else if(stat_b.pmin != st.pmin || stat_b.pmax != st.pmax || stat_b.nmin != st.nmin || stat_b.nmax != st.nmax ||
			stat_b.kmin != st.kmin || stat_b.kmax != st.kmax || stat_b.cw != st.cw){
			fprintf(stderr,"Invalid checkpoint file %s !!!\n",STATE_FILENAME_B);
			printf("Invalid checkpoint file %s !!!\n",STATE_FILENAME_B);
			good_state_b = false;
		}
		else{
			uint64_t state_sum = stat_b.pmin+stat_b.pmax+stat_b.p+stat_b.checksum+stat_b.primecount+stat_b.factorcount+
						stat_b.last_trickle+stat_b.nmin+stat_b.nmax+stat_b.kmin+stat_b.kmax+stat_b.cw;

			if(state_sum != stat_b.state_sum){
				fprintf(stderr,"Checksum error in %s !!!\n",STATE_FILENAME_B);
				printf("Checksum error in %s !!!\n",STATE_FILENAME_B);
				good_state_b = false;
			}
		}
		fclose(in);
	}

        // If both state files are OK, check which is the most recent
	if (good_state_a && good_state_b)
	{
		if (stat_a.p > stat_b.p)
			good_state_b = false;
		else
			good_state_a = false;
	}

        // Use data from the most recent state file
	if (good_state_a && !good_state_b)
	{
		memcpy(&st, &stat_a, sizeof(workStatus));
		sd.write_state_a_next = false;
		if(boinc_is_standalone()){
			printf("Resuming from checkpoint in %s\n",STATE_FILENAME_A);
		}
		return 1;
	}
        if (good_state_b && !good_state_a)
        {
		memcpy(&st, &stat_b, sizeof(workStatus));
		sd.write_state_a_next = true;
		if(boinc_is_standalone()){
			printf("Resuming from checkpoint in %s\n",STATE_FILENAME_B);
		}
		return 1;
        }

	// If we got here, neither state file was good
	return 0;
}

void checkpoint( workStatus & st, searchData & sd ){
	handle_trickle_up( st );
	write_state( st, sd );
	if(boinc_is_standalone()){
		printf("Checkpoint, current p: %" PRIu64 "\n", st.p);
	}
	boinc_checkpoint_completed();
}


// sleep CPU thread while waiting on the specified event to complete in the command queue
void waitOnEvent(sclHard hardware, cl_event event)
{
    cl_int err;
    cl_int info;

#ifndef _WIN32
    struct timespec sleep_time;
    sleep_time.tv_sec  = 0;
    sleep_time.tv_nsec = 1000000; // 1ms
#endif

    boinc_begin_critical_section();

    err = clFlush(hardware.queue);
    if (err != CL_SUCCESS) {
        printf("ERROR: clFlush\n");
        fprintf(stderr, "ERROR: clFlush\n");
        sclPrintErrorFlags(err);
    }

    while (true) {

#ifdef _WIN32
	Sleep(1);
#else
        nanosleep(&sleep_time, NULL);
#endif

        err = clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS,
                             sizeof(cl_int), &info, NULL);
        if (err != CL_SUCCESS) {
            printf("ERROR: clGetEventInfo\n");
            fprintf(stderr, "ERROR: clGetEventInfo\n");
            sclPrintErrorFlags(err);
        }

        if (info == CL_COMPLETE) {
            err = clReleaseEvent(event);
            if (err != CL_SUCCESS) {
                printf("ERROR: clReleaseEvent\n");
                fprintf(stderr, "ERROR: clReleaseEvent\n");
                sclPrintErrorFlags(err);
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
#ifndef _WIN32
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

int factorcompare(const void *a, const void *b) {
  	factor *factA = (factor *)a;
	factor *factB = (factor *)b;
	if(factB->p < factA->p){
		return 1;
	}
	else if(factB->p == factA->p){
		if(factB->n < factA->n){
			return 1;
		}
	}
	return -1;
}

void getResults( progData & pd, workStatus & st, searchData & sd, sclHard hardware, cl_ulong *h_checksum, cl_uint *h_primecount ){

	if(boinc_is_standalone()){
		printf("\r                                                                                \r");
		fflush(stdout);
	}

	// copy checksum and total prime count to host memory, non-blocking
	sclReadNB(hardware, sd.numgroups*sizeof(cl_ulong), pd.d_checksum, h_checksum);
	// copy prime count to host memory, blocking
	sclRead(hardware, 4*sizeof(cl_uint), pd.d_primecount, h_primecount);

	// index 0 is the gpu's total prime count
	st.primecount += h_checksum[0];

	// sum block checksums
	for(uint32_t i=1; i<sd.numgroups; ++i){
		st.checksum += h_checksum[i];
	}

	// largest kernel prime count.  used to check array bounds
	if(h_primecount[1] > sd.psize){
		fprintf(stderr,"error: gpu prime array overflow\n");
		printf("error: gpu prime array overflow\n");
		exit(EXIT_FAILURE);
	}

	// flag set by gpu if there is an internal checksum error
	if(h_primecount[3]){
		fprintf(stderr,"error: gpu checksum failure\n");
		printf("error: gpu checksum failure\n");
		exit(EXIT_FAILURE);
	}

	uint32_t numfactors = h_primecount[2];
	if(numfactors){

		if(numfactors > sd.numresults){
			fprintf(stderr,"Error: number of results (%u) overflowed array.\n", numfactors);
			printf("Error: number of results (%u) overflowed array.\n", numfactors);
			exit(EXIT_FAILURE);
		}

		factor *h_factor = (factor *)malloc(numfactors * sizeof(factor));
		if( h_factor == NULL ){
			fprintf(stderr,"malloc error: h_factor\n");
			printf("malloc error: h_factor\n");
			exit(EXIT_FAILURE);
		}

		// copy factors to host memory
		// blocking read
		sclRead(hardware, numfactors * sizeof(factor), pd.d_factor, h_factor);

		// sort factors by prime size if needed
		if(numfactors > 1){
			qsort(h_factor, numfactors, sizeof(factor), factorcompare);
		}

		FILE * resfile = NULL;

		if(boinc_is_standalone()){
			printf("Verifying %u factors on CPU...\n", numfactors);
		}

		for(uint32_t i=0; i<numfactors; ++i){
			uint64_t fp = h_factor[i].p;
			uint32_t fn = h_factor[i].n;
			uint32_t fk = (h_factor[i].k < 0) ? -h_factor[i].k : h_factor[i].k;
			int32_t fc = (h_factor[i].k < 0) ? -1 : 1;

			// from ppsieve, not needed?
			if(!st.cw){
				if( (fk&1)==0 ) continue;	// k is even
			}

			// check for a small prime factor of the number
			// is this factor useful?
			uint64_t smf = try_all_factors(fk, fn, fc); 
			if(!smf){
				// no factors under 2^15	
				// check the large factor actually divides the number
				int32_t vres = verify_factor(fp,fk,fn,fc); 
				if( !vres ){
					fprintf(stderr,"CPU factor verification failed!  %" PRIu64 " is not a factor of %u*2^%u%+d\n", fp, fk, fn, fc);
					printf("CPU factor verification failed!  %" PRIu64 " is not a factor of %u*2^%u%+d\n", fp, fk, fn, fc);
					exit(EXIT_FAILURE);
				}
				else if( vres == -1 ){
					// Unlikely, factor is a 2-prp
//					++prpcount;
				}
				else{
					if(resfile == NULL){
						resfile = my_fopen(RESULTS_FILENAME,"a");
						if( resfile == NULL ){
							fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
							printf("Cannot open %s !!!\n",RESULTS_FILENAME);
							exit(EXIT_FAILURE);
						}
					}
					
					if( fprintf( resfile, "%" PRIu64 " | %u*2^%u%+d\n",fp,fk,fn,fc) < 0 ){
						fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
						printf("Cannot write to %s !!!\n",RESULTS_FILENAME);
						exit(EXIT_FAILURE);
					}						
					// add the factor to checksum
					st.checksum += fk + fn + fc;
					++st.factorcount;
				}
			}
			else{
//				printf("factor: %" PRIu64 " | %u*2^%u%+d has a small prime factor %" PRIu64 "\n",fp,fk,fn,fc,smf);
			}

		}

		if(resfile != NULL){
			fclose(resfile);
		}

		free(h_factor);
	}
	checkpoint(st,sd);
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


void setupSearch(workStatus & st, searchData & sd){

	// increase nmax to check for factors equal to nmax
	++st.nmax;

	st.p = st.pmin;

	if(st.pmin == 0 || st.pmax == 0){
		printf("-p and -P arguments are required\n");
		fprintf(stderr, "-p and -P arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if(st.nmin == 0 || st.nmax == 0){
		printf("-n and -N arguments are required\n");
		fprintf(stderr, "-n and -N arguments are required\n");
		exit(EXIT_FAILURE);
	}

	if (st.nmin > st.nmax){
		printf("nmin <= nmax is required\n");
		fprintf(stderr, "nmin <= nmax is required\n");
		exit(EXIT_FAILURE);
	}

	if(st.cw){

		if(st.nmax >= st.pmin){
			printf("nmax < pmin is required\n");
			fprintf(stderr, "nmax < pmin is required\n");
			exit(EXIT_FAILURE);
		}

		st.kmax = st.nmax;
		st.kmin = st.nmin;
	}
	else{

		if(st.kmax == 0){
			printf("-K argument is required\n");
			fprintf(stderr, "-K argument is required\n");
			exit(EXIT_FAILURE);
		}

		if(st.kmin > st.kmax){
			printf("kmin <= kmax is required\n");
			fprintf(stderr, "kmin <= kmax is required\n");
			exit(EXIT_FAILURE);
		}

		if(st.kmax >= st.pmin){
			printf("kmax < pmin is required\n");
			fprintf(stderr, "kmax < pmin is required\n");
			exit(EXIT_FAILURE);
		}

		uint32_t b0 = 0, b1 = 0;
		b0 = st.kmin/sd.kstep;
		b1 = st.kmax/sd.kstep;
		st.kmin = b0*sd.kstep+sd.koffset;
		st.kmax = b1*sd.kstep+sd.koffset;
	}


	for (sd.nstep = 1; ( (uint64_t)(st.kmax) << sd.nstep ) < st.pmin; sd.nstep++);

	if((((uint64_t)1) << (64-sd.nstep)) > st.pmin) {

		uint64_t pmin_1 = (((uint64_t)1) << (64-sd.nstep));

		printf("Error: pmin is not large enough (or nmax is close to nmin).\n");
		fprintf(stderr, "Error: pmin is not large enough (or nmax is close to nmin).\n");

		st.pmin = st.kmax + 1;
		for (sd.nstep = 1; ( (uint64_t)(st.kmax) << sd.nstep ) < st.pmin; sd.nstep++);

		while((((uint64_t)1) << (64-sd.nstep)) > st.pmin) {
			st.pmin *= 2;
			sd.nstep++;
		}
		if(pmin_1 < st.pmin){
			st.pmin = pmin_1;
		}

		printf("This program will work by the time pmin == %" PRIu64 ".\n", st.pmin);
		fprintf(stderr, "This program will work by the time pmin == %" PRIu64 ".\n", st.pmin);

		exit(EXIT_FAILURE);
	}

	if (sd.nstep > (st.nmax-st.nmin+1))
		sd.nstep = (st.nmax-st.nmin+1);

	// For TPS, decrease the ld_nstep by one to allow overlap, checking both + and -
	sd.nstep--;

	// Use the 32-step algorithm where useful.
	if(sd.nstep >= 32 && (((uint64_t)1) << 32) <= st.pmin) {
		sd.nstep = 32;
	}

// for testing
// sd.nstep = 20;

	// search twin, decrease by one to allow overlap, checking both + and -
	st.nmin--;

	sd.bbits = lg2(st.nmin);

	if(sd.bbits < 6) {
		printf("Error: nmin too small at %d (must be at least 65).\n", st.nmin+1);
		fprintf(stderr, "Error: nmin too small at %d (must be at least 65).\n", st.nmin+1);
		exit(EXIT_FAILURE);
	}

	// r = 2^-i * 2^64 (mod N), something that can be done in a uint64_t!
	// If i is large (and it should be at least >= 32), there's a very good chance no mod is needed!
	sd.r0 = ((uint64_t)1) << (64-(st.nmin >> (sd.bbits-5)));

	sd.bbits = sd.bbits-6;

	sd.mont_nstep = 64-sd.nstep;

	// data for checksum
	uint32_t maxn;

	maxn = ( (st.nmax-st.nmin) / sd.nstep) * sd.nstep;
	maxn += st.nmin;

	if( maxn < st.nmax ){
		maxn += sd.nstep;
	}

	int bbits1 = lg2(maxn) - 5;
	sd.r1 = ((uint64_t)1) << (64-(maxn >> bbits1));
	--bbits1;
	sd.bbits1 = bbits1;
	sd.lastN = maxn;

	// increase result buffer at low P range
	// it's still possible to overflow this with a fast GPU and large search range
	if(st.pmin < 0xFFFFFFFF){
		sd.numresults = 30000000;
	}
	else{
		sd.numresults = 10000000;
	}
}


void profileGPU(progData & pd, workStatus & st, searchData & sd, sclHard hardware, int debuginfo ){

	// calculate approximate chunk size based on gpu's compute units
	cl_int err = 0;

	uint64_t calc_range = sd.computeunits * 2000000;

	// limit kernel global size
	if(calc_range > 4294900000){
		calc_range = 4294900000;
	}

	uint64_t estimated = calc_range;

	uint64_t prof_start = st.p;

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

	// target is 5ms for prime generator kernel
	double prof_multi = 5.0 / kernel_ms;

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

	sd.range = calc_range;
	sd.psize = mem_size;

	// free temporary array
	sclReleaseMemObject(d_profileprime);

	// N's to search each time a kernel is run
	if(sd.compute){
		sd.kernel_nstep = sd.nstep * 3750;
	}
	else{
		sd.kernel_nstep = sd.nstep * 750;
	}

}


void cl_sieve( sclHard hardware, searchData & sd, workStatus & st ){

	progData pd;
	bool profile = true;
	bool debuginfo = false;
	cl_int err = 0;

	sieve_small_primes();

	// setup kernel parameters
	setupSearch(st,sd);

	fprintf(stderr, "Starting sieve at p: %" PRIu64 " n: %u k: %u\nStopping sieve at P: %" PRIu64 " N: %u K: %u\n", st.pmin, st.nmin+1, st.kmin, st.pmax, st.nmax, st.kmax);
	if(boinc_is_standalone()){
		printf("Starting sieve at p: %" PRIu64 " n: %u k: %u\nStopping sieve at P: %" PRIu64 " N: %u K: %u\n", st.pmin, st.nmin+1, st.kmin, st.pmax, st.nmax, st.kmax);
	}

	// device arrays
	pd.d_primecount = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, 4*sizeof(cl_uint), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: d_primecount array.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
        pd.d_factor = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numresults*sizeof(factor), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	// host
	cl_uint * h_primecount = (cl_uint*)malloc(4 * sizeof(cl_uint));
	if( h_primecount == NULL ){
		fprintf(stderr,"malloc error: h_primecount\n");
		exit(EXIT_FAILURE);
	}

	// bake constants into kernel source
	char cldef[256];
	snprintf(cldef, sizeof(cldef), "-DNSTEP=%u -DMONT_NSTEP=%u -DNMAX=%u -DKMIN=%u -DKMAX=%u -DSM_MONT_NSTEP=%u",
		sd.nstep, sd.mont_nstep, st.nmax, st.kmin, st.kmax, sd.mont_nstep-32);
//	printf("%s\n", cldef);
	if(st.cw){
		if(sd.nstep == 32){
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecw32",hardware,cldef);
		}
		else if(sd.nstep < 32){
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecwsm",hardware,cldef);
		}
		else{
			pd.sieve = sclGetCLSoftware(sievecw_cl,"sievecw",hardware,cldef);
		}
	}
	else{
		if(sd.nstep == 32){
			pd.sieve = sclGetCLSoftware(sieve_cl,"sieve32",hardware,cldef);
		}
		else if(sd.nstep < 32){
			pd.sieve = sclGetCLSoftware(sieve_cl,"sievesm",hardware,cldef);
		}
		else{
			pd.sieve = sclGetCLSoftware(sieve_cl,"sieve",hardware,cldef);
		}
	}

        pd.clearn = sclGetCLSoftware(clearn_cl,"clearn",hardware,NULL);
        pd.clearresult = sclGetCLSoftware(clearresult_cl,"clearresult",hardware,NULL);
        pd.check = sclGetCLSoftware(check_cl,"check",hardware,NULL);
        pd.getsegprimes = sclGetCLSoftware(getsegprimes_cl,"getsegprimes",hardware,NULL);

	// bake constants into kernel source
	snprintf(cldef, sizeof(cldef), "-DNMIN=%u -DLASTN=%u -DBBITS=%d -DBBITSL=%d -DRS=((ulong)%" PRIu64 ") -DRSL=((ulong)%" PRIu64 ")",
    					st.nmin, sd.lastN, sd.bbits, sd.bbits1, sd.r0, sd.r1);
//	printf("%s\n", cldef);
        pd.setup = sclGetCLSoftware(setup_cl,"setup",hardware,cldef);

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
		if( read_state( st, sd ) ){
			if(boinc_is_standalone()){
				printf("Resuming search from checkpoint. Current p: %" PRIu64 "\n", st.p);
			}
			fprintf(stderr,"Resuming search from checkpoint. Current p: %" PRIu64 "\n", st.p);

			//trying to resume a finished workunit
			if( st.p == st.pmax ){
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
			st.last_trickle = (uint64_t)time(NULL);
		}
	}

	// kernel used in profileGPU, setup arg
	sclSetKernelArg(pd.clearn, 0, sizeof(cl_mem), &pd.d_primecount);
	sclSetGlobalSize( pd.clearn, 64 );

	profileGPU(pd,st,sd,hardware,debuginfo);

	// number of gpu workgroups, used to size the checksum array on gpu
	sd.numgroups = (sd.psize / pd.check.local_size[0]) + 2;

	sclSetGlobalSize( pd.getsegprimes, (sd.range/60)+1 );
	sclSetGlobalSize( pd.setup, sd.psize );
	sclSetGlobalSize( pd.sieve, sd.psize );
	sclSetGlobalSize( pd.check, sd.psize );
	sclSetGlobalSize( pd.clearresult, sd.numgroups );

	// allocate gpu P, Ps, K, lastK arrays
	pd.d_primes = clCreateBuffer(hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong), NULL, &err);
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_Ps = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_K = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	pd.d_lK = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.psize*sizeof(cl_ulong), NULL, &err );
	if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}

        pd.d_checksum = clCreateBuffer( hardware.context, CL_MEM_READ_WRITE, sd.numgroups*sizeof(cl_ulong), NULL, &err );
        if ( err != CL_SUCCESS ) {
		fprintf(stderr, "ERROR: clCreateBuffer failure: d_checksum array.\n");
                printf( "ERROR: clCreateBuffer failure.\n" );
		exit(EXIT_FAILURE);
	}
	// host
	cl_ulong *h_checksum = (cl_ulong*)malloc(sd.numgroups * sizeof(cl_ulong));
	if( h_primecount == NULL ){
		fprintf(stderr,"malloc error: h_primecount\n");
		exit(EXIT_FAILURE);
	}

	// set static kernel args
	int ai=0;
	sclSetKernelArg(pd.clearresult, ai++, sizeof(cl_mem), &pd.d_checksum);
	sclSetKernelArg(pd.clearresult, ai++, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.clearresult, ai++, sizeof(uint32_t), &sd.numgroups);
	////////////////////////
	sclSetKernelArg(pd.getsegprimes, 3, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.getsegprimes, 4, sizeof(cl_mem), &pd.d_primecount);
	////////////////////////
	ai=0;
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_Ps);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_lK);
	sclSetKernelArg(pd.setup, ai++, sizeof(cl_mem), &pd.d_primecount);
	////////////////////////
	ai=0;
	sclSetKernelArg(pd.sieve, ai++, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.sieve, ai++, sizeof(cl_mem), &pd.d_Ps);
	sclSetKernelArg(pd.sieve, ai++, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.sieve, ai++, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.sieve, ai++, sizeof(cl_mem), &pd.d_factor);
	////////////////////////
	ai=0;
	sclSetKernelArg(pd.check, ai++, sizeof(cl_mem), &pd.d_K);
	sclSetKernelArg(pd.check, ai++, sizeof(cl_mem), &pd.d_lK);
	sclSetKernelArg(pd.check, ai++, sizeof(cl_mem), &pd.d_primecount);
	sclSetKernelArg(pd.check, ai++, sizeof(cl_mem), &pd.d_primes);
	sclSetKernelArg(pd.check, ai++, sizeof(cl_mem), &pd.d_checksum);
	sclSetKernelArg(pd.check, ai++, sizeof(uint32_t), &sd.numgroups);
	////////////////////////

	fprintf(stderr,"Starting sieve...\n");
	if(boinc_is_standalone()){
		printf("Starting sieve...\n");
	}

	time_t boinc_last, ckpt_last, time_curr;
	time(&boinc_last);
	time(&ckpt_last);

	// clear results, checksum, total prime counts
	sclEnqueueKernel(hardware, pd.clearresult);

	time_t totals, totalf;
	if(boinc_is_standalone()){
		time(&totals);
	}

	const double irsize = 1.0 / (double)(st.pmax-st.pmin);
	time_t start_time = time(NULL);
	double smooth_rate = 0;
	uint64_t run_start_p = st.p;

	// main search loop
	for(uint64_t stop; st.p < st.pmax; st.p = stop){

		// clear prime count
		sclEnqueueKernel(hardware, pd.clearn);

		stop = st.p + sd.range;
		if(stop > st.pmax) stop = st.pmax;

		time(&time_curr);
		if( ((int)time_curr - (int)boinc_last) > 1 ){
			// update BOINC fraction done every 2 sec
    			double fd = (double)(st.p-st.pmin)*irsize;
			boinc_fraction_done(fd);
			if(boinc_is_standalone()){
				print_progress(st, &smooth_rate, start_time, run_start_p);
			}
			boinc_last = time_curr;
			int elapsed = (int)time_curr - (int)ckpt_last;
			if(elapsed > 60){
				sleepCPU(hardware);
				boinc_begin_critical_section();
				getResults(pd, st, sd, hardware, h_checksum, h_primecount);
				boinc_end_critical_section();
				ckpt_last = time_curr;
				// clear result arrays
				sclEnqueueKernel(hardware, pd.clearresult);	
			}
		}

		// get primes
		int32_t wheelidx;
		uint64_t kernel_start = st.p;
		findWheelOffset(kernel_start, wheelidx);

		sclSetKernelArg(pd.getsegprimes, 0, sizeof(uint64_t), &kernel_start);
		sclSetKernelArg(pd.getsegprimes, 1, sizeof(uint64_t), &stop);
		sclSetKernelArg(pd.getsegprimes, 2, sizeof(int32_t), &wheelidx);
		cl_event launchEvent = sclEnqueueKernelEvent(hardware, pd.getsegprimes);

		// setup Ps, K kernel
		sclEnqueueKernel(hardware, pd.setup);

		uint32_t nstart = st.nmin;
		uint32_t nend;

		// profile gpu sieve kernel time once, at program start.  adjust work size to target kernel runtime.
		if(profile){
			nend = nstart + sd.kernel_nstep;
			if(nend > st.nmax) nend = st.nmax;
			sclSetKernelArg(pd.sieve, 5, sizeof(uint32_t), &nstart);
			sclSetKernelArg(pd.sieve, 6, sizeof(uint32_t), &nend);
			double kernel_ms = ProfilesclEnqueueKernel(hardware, pd.sieve);
			double multi = (sd.compute)?(50.0 / kernel_ms):(10.0 / kernel_ms);	// target kernel time 50ms or 10ms
			uint32_t new_knstep = (uint32_t)((double)sd.kernel_nstep * multi);
			// make sure it's a multiple of nstep
			new_knstep = (new_knstep / sd.nstep) * sd.nstep;
			if(debuginfo) printf("old kns %u, new kns %u\n",sd.kernel_nstep,new_knstep);
			sd.kernel_nstep = new_knstep;
			fprintf(stderr,"c:%u u:%u r:%u p:%u ns:%u kns:%u\n", (uint32_t)sd.compute, sd.computeunits, sd.range, sd.psize, sd.nstep, sd.kernel_nstep);
			if(boinc_is_standalone()){
				printf("c:%u u:%u r:%u p:%u ns:%u kns:%u\n", (uint32_t)sd.compute, sd.computeunits, sd.range, sd.psize, sd.nstep, sd.kernel_nstep);
			}
			profile = false;
			nstart = nend;
		}

		// sieve kernel, loop to nmax
		while(nstart < st.nmax){
			nend = nstart + sd.kernel_nstep;
			if(nend > st.nmax) nend = st.nmax;
			sclSetKernelArg(pd.sieve, 5, sizeof(uint32_t), &nstart);
			sclSetKernelArg(pd.sieve, 6, sizeof(uint32_t), &nend);
			sclEnqueueKernel(hardware, pd.sieve);
//			float kernel_ms = ProfilesclEnqueueKernel(hardware, pd.sieve);
//			printf("sieve kernel time %0.2fms\n",kernel_ms);
			nstart = nend;
		}

		// validate checksum kernel
		sclEnqueueKernel(hardware, pd.check);

		// limit cl queue depth and sleep cpu
		waitOnEvent(hardware, launchEvent);

	}


	// final checkpoint
	sleepCPU(hardware);
	boinc_begin_critical_section();
	st.p = st.pmax;
	boinc_fraction_done(1.0);
	getResults(pd, st, sd, hardware, h_checksum, h_primecount);
	
	// print checksum
	FILE * resfile = my_fopen(RESULTS_FILENAME,"a");

	if(resfile == NULL){
		fprintf(stderr,"Cannot open %s !!!\n",RESULTS_FILENAME);
		exit(EXIT_FAILURE);
	}

	if(st.factorcount){
		if( fprintf( resfile, "%016" PRIX64 "\n", st.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}
	else{
		if( fprintf( resfile, "no factors\n%016" PRIX64 "\n", st.checksum ) < 0 ){
			fprintf(stderr,"Cannot write to %s !!!\n",RESULTS_FILENAME);
			exit(EXIT_FAILURE);
		}
	}

	fclose(resfile);
	
	boinc_end_critical_section();

	fprintf(stderr,"Sieve complete.\nfactors %" PRIu64 ", prime count %" PRIu64 "\n", st.factorcount, st.primecount);

	if(boinc_is_standalone()){
		time(&totalf);
		printf("Sieve finished in %d sec.\n", (int)totalf - (int)totals);
		printf("factors %" PRIu64 ", prime count %" PRIu64 ", checksum %016" PRIX64 "\n", st.factorcount, st.primecount, st.checksum);
	}

	// cleanup
	free(h_primecount);
	free(h_checksum);
	cleanup(pd);
	small_primes_free();
}


void run_test( sclHard hardware, searchData & sd, workStatus & st){

	int goodtest = 0;

	printf("Beginning self test of 6 ranges.\n");

//	-p 25636026e6 -P 25636030e6 -n 10000000 -N 25000000 -c		nstep 19
	st.pmin = 25636026000000;
	st.pmax = 25636030000000;
	st.nmin = 10000000;
	st.nmax = 25000000;
	st.kmin = 0;
	st.kmax = 0;
	st.cw = true;
	cl_sieve( hardware, sd, st );
	if( st.factorcount == 2 && st.primecount == 129869 && st.checksum == 0x4544591DC69ACD83 ){
		printf("CW test case 1 passed.\n\n");
		fprintf(stderr,"CW test case 1 passed.\n");
		++goodtest;
	}
	else{
		printf("CW test case 1 failed.\n\n");
		fprintf(stderr,"CW test case 1 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;

//	-p 556439300e6 -P 556439440e6 -n 100 -N 100000 -c		nstep 32
	st.pmin = 556439300000000;
	st.pmax = 556439440000000;
	st.nmin = 100;
	st.nmax = 100000;
	st.kmin = 0;
	st.kmax = 0;
	st.cw = true;
	cl_sieve( hardware, sd, st);
	if( st.factorcount == 1 && st.primecount == 4123452 && st.checksum == 0x8FEC30979896A3C0 ){
		printf("CW test case 2 passed.\n\n");
		fprintf(stderr,"CW test case 2 passed.\n");
		++goodtest;
	}
	else{
		printf("CW test case 2 failed.\n\n");
		fprintf(stderr,"CW test case 2 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;


//	-p838338347800e6 -P838338347820e6 -k5 -K9999 -n6000000 -N9000000	nstep 32
	st.pmin = 838338347800000000;
	st.pmax = 838338347820000000;
	st.nmin = 6000000;
	st.nmax = 9000000;
	st.kmin = 5;
	st.kmax = 9999;
	st.cw = false;
	cl_sieve( hardware, sd, st );
	if( st.factorcount == 1 && st.primecount == 484024 && st.checksum == 0xA7DC855BCB311759 ){
		printf("test case 3 passed.\n\n");
		fprintf(stderr,"test case 3 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 3 failed.\n\n");
		fprintf(stderr,"test case 3 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;

//	-p42070000e6 -P42070050e6 -k 1201 -K 9999 -n 100 -N 2000000		nstep 31
	st.pmin = 42070000000000;
	st.pmax = 42070050000000;
	st.nmin = 100;
	st.nmax = 2000000;
	st.kmin = 1201;
	st.kmax = 9999;
	st.cw = false;
	cl_sieve( hardware, sd, st );
	if( st.factorcount == 70 && st.primecount == 1592285 && st.checksum == 0x727796B2D3677937 ){
		printf("test case 4 passed.\n\n");
		fprintf(stderr,"test case 4 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 4 failed.\n\n");
		fprintf(stderr,"test case 4 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;


//	-k 5 -K 9999 -n 65 -N 3000000 -p 47772822600000 -P 47773822700000	nstep 32
	st.pmin = 47772822600000;
	st.pmax = 47773822700000;
	st.nmin = 65;
	st.nmax = 3000000;
	st.kmin = 5;
	st.kmax = 9999;
	st.cw = false;
	cl_sieve( hardware, sd, st );
	if( st.factorcount == 2255 && st.primecount == 31755968 && st.checksum == 0x5DE2E7801F431850 ){
		printf("test case 5 passed.\n\n");
		fprintf(stderr,"test case 5 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 5 failed.\n\n");
		fprintf(stderr,"test case 5 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;

//	-k 5042 -K 99999 -n 15124185 -N 20001342 -p 11e11 -P 110001e7	nstep 23
	st.pmin = 1100000000000;
	st.pmax = 1100010000000;
	st.nmin = 15124185;
	st.nmax = 20001342;
	st.kmin = 5042;
	st.kmax = 99999;
	st.cw = false;
	cl_sieve( hardware, sd, st );
	if( st.factorcount == 16286 && st.primecount == 361226 && st.checksum == 0x08466535C212C86B ){
		printf("test case 6 passed.\n\n");
		fprintf(stderr,"test case 6 passed.\n");
		++goodtest;
	}
	else{
		printf("test case 6 failed.\n\n");
		fprintf(stderr,"test case 6 failed.\n");
	}
	st.checksum = 0;
	st.primecount = 0;
	st.factorcount = 0;


	if(goodtest == 6){
		printf("All test cases completed successfully!\n");
		fprintf(stderr, "All test cases completed successfully!\n");
	}
	else{
		printf("Self test FAILED!\n");
		fprintf(stderr, "Self test FAILED!\n");
	}

}


