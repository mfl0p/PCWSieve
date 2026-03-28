/*
	Montgomery functions

	Bryan Little, Mar 2026
*/

#include "mont.h"

uint64_t invert(uint64_t p){
	uint64_t p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

uint64_t m_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t q){
	unsigned __int128 res;
	res  = (unsigned __int128)a * b;
	uint64_t ab0 = (uint64_t)res;
	uint64_t ab1 = res >> 64;
	uint64_t m = ab0 * q;
	res = (unsigned __int128)m * p;
	uint64_t mp = res >> 64;
	uint64_t r = ab1 - mp;
	return ( ab1 < mp ) ? r + p : r;
}

uint64_t add(uint64_t a, uint64_t b, uint64_t p){
	uint64_t r;
	uint64_t c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}
