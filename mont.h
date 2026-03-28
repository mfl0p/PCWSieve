/*
	Montgomery functions

	Bryan Little, Mar 2026
*/

#ifndef MONT_H
#define MONT_H

#include <cinttypes>

uint64_t invert(uint64_t p);

uint64_t m_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t q);

uint64_t add(uint64_t a, uint64_t b, uint64_t p);

#endif
