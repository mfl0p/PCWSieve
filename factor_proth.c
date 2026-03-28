/* 
	factor_proth.c

	Bryan Little, Mar 2026

	montgomery arithmetic by Yves Gallot

	Factor a Proth number with small primes, and see if it breaks.
	To be used to test whether potential larger factors are useful.
*/

#include "mont.h"
#include "factor_proth.h"
#include "primesieve.h"

uint64_t* small_primes;
int32_t spcount;

// Returns a list of small primes, 3-32767 inclusive.
void sieve_small_primes() {
	size_t size;
	small_primes = (uint64_t*) primesieve_generate_primes(3, 32767, &size, UINT64_PRIMES);
	spcount = (int32_t)size;
}

void small_primes_free(){
	primesieve_free(small_primes);
}

int32_t try_factor(uint32_t k, uint32_t n, int32_t c, uint64_t p) {
	uint64_t q = invert(p);
	uint64_t one = (-p) % p;
	uint64_t two = add(one, one, p);
	uint32_t exp = n;
	uint32_t curBit = 0x80000000;
	curBit >>= ( __builtin_clz(exp) + 1 );

	uint64_t a = two;

	// a = 2^n mod P
	while( curBit ){
		a = m_mul(a,a,p,q);
		if(exp & curBit){
			a = add(a,a,p);
		}
		curBit >>= 1;
	}
 
	// a = k*2^n mod P, not in montgomery form
	a = m_mul(a,k,p,q);

	if(a == 1 && c == -1){
//		printf("%" PRIu64 " is a factor of %u*2^%u-1\n",p,k,n);
		return 1;
	}
	else if(a == p-1 && c == 1){
//		printf("%" PRIu64 " is a factor of %u*2^%u+1\n",p,k,n);
		return 1;
	}
	return 0;
}

uint64_t try_all_factors(uint32_t k, uint32_t n, int32_t c) {
	for(int i=0; i < spcount; ++i){
		if(try_factor(k, n, c, small_primes[i])){
			return small_primes[i];
		}
	}
	return 0;
}


