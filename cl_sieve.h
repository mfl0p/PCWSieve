
// cl_sieve.h

typedef struct {
	uint64_t pmin, pmax, checksum, r0, r1, lastN, p, workunit, last_trickle, primecount, factorcount;
	uint32_t kmin, kmax, nmin, nmax, nstep, mont_nstep, kernel_nstep, numresults, kstep, koffset;
	int32_t bbits, bbits1, computeunits;
	bool write_state_a_next, cw, test, compute;
}searchData;

void cl_sieve( sclHard hardware, searchData & sd );

void run_test( sclHard hardware, searchData & sd );
