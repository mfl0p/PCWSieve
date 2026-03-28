
// cl_sieve.h

typedef struct {
	cl_ulong p;
	cl_uint n;
	cl_int k;
} factor;

typedef struct {
	uint64_t pmin, pmax, p, checksum, primecount, factorcount, last_trickle, state_sum;
	uint32_t nmin, nmax, kmin, kmax;
	bool cw;
} workStatus;

typedef struct {
	uint64_t r0, r1;
	uint32_t nstep, mont_nstep, kernel_nstep, numresults, kstep, koffset, numgroups, psize, range, lastN;
	int32_t bbits, bbits1, computeunits;
	bool write_state_a_next, test, compute;
} searchData;

typedef struct {
	cl_mem d_factor = NULL;
	cl_mem d_checksum = NULL;

	cl_mem d_primes = NULL;
	cl_mem d_primecount = NULL;

	cl_mem d_Ps = NULL;
	cl_mem d_K = NULL;
	cl_mem d_lK = NULL;

	sclSoft sieve, clearn, clearresult, setup, check, getsegprimes;

}progData;

void cl_sieve( sclHard hardware, searchData & sd, workStatus & st );

void run_test( sclHard hardware, searchData & sd, workStatus & st );
