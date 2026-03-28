/* 
	factor_proth.h 
	Factor a Proth number with small primes, and see if it breaks.
	To be used to test whether potential larger factors are useful.
*/

#ifndef FACTOR_PROTH_H
#define FACTOR_PROTH_H

void sieve_small_primes();

void small_primes_free();

uint64_t try_all_factors(uint32_t k, uint32_t n, int32_t c);

#endif
