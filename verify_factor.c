/* 
	verify_factor.c

	Bryan Little, Mar 2026

	montgomery arithmetic by Yves Gallot
*/

#include "mont.h"
#include "verify_factor.h"

/* Used in the prime validator
   Returns 0 only if p is composite.
   Otherwise p is a strong probable prime to base a.
 */
bool strong_prp(uint32_t base, uint64_t p, uint64_t q, uint64_t one, uint64_t pmo, uint64_t r2, int t, uint64_t exp, uint64_t curBit){

	/* If p is prime and p = d*2^t+1, where d is odd, then either
		1.  a^d = 1 (mod p), or
		2.  a^(d*2^s) = -1 (mod p) for some s in 0 <= s < t    */

	uint64_t a = m_mul(base,r2,p,q);  // convert base to montgomery form
	const uint64_t mbase = a;

  	/* r <-- a^d mod p, assuming d odd */
	while( curBit )
	{
		a = m_mul(a,a,p,q);

		if(exp & curBit){
			a = m_mul(a,mbase,p,q);
		}

		curBit >>= 1;
	}

	/* Clause 1. and s = 0 case for clause 2. */
	if (a == one || a == pmo){
		return true;
	}

	/* 0 < s < t cases for clause 2. */
	for (int s = 1; s < t; ++s){

		a = m_mul(a,a,p,q);

		if(a == pmo){
	    		return true;
		}
	}


	return false;
}


// prime if the number passes prp test to 7 bases.  good to 2^64
// this is very fast
bool isPrime(uint64_t p, uint64_t q, uint64_t one, uint64_t pmo, uint64_t two, uint64_t r2){

	const uint32_t bases[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

	if (p % 2==0)
		return false;

	int t = __builtin_ctzll( (p-1) );
	uint64_t exp = p >> t;
	uint64_t curBit = 0x8000000000000000;
	curBit >>= ( __builtin_clzll(exp) + 1 );

	for (int i = 0; i < 7; ++i){

		uint32_t base = bases[i];

		// needed for composite bases
		if (base >= p){
			base %= p;
			if (base == 0)
				continue;
		}

		if (!strong_prp(base, p, q, one, pmo, r2, t, exp, curBit))
			return false;
	}

	return true;
}

int verify_factor(uint64_t p, uint32_t k, uint32_t n, int32_t c){

	uint64_t q = invert(p);
	uint64_t one = (-p) % p;
	uint64_t pmo = p - one;	
	uint64_t two = add(one, one, p);
	uint64_t t = add(two, two, p);
	for (int i = 0; i < 5; ++i)
		t = m_mul(t, t, p, q);	// 4^{2^5} = 2^64
	uint64_t r2 = t;

	if(!isPrime(p, q, one, pmo, two, r2)){
		return -1;
	}

	uint32_t exp = n;
	uint32_t curBit = 0x80000000;
	curBit >>= ( __builtin_clz(exp) + 1 );

	uint64_t a = two;

	// a = 2^n mod P
	while( curBit )	{
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



