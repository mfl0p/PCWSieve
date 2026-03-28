/*

	setup kernel

	Bryan Little 2/12/2023
	Ken Brazier August 2010

	setup data for sieve and check kernels

*/


ulong mulmod_REDC (const ulong a, const ulong b, const ulong N, const ulong Ns)
{
        ulong rax, rcx;

#ifdef __NV_CL_C_VERSION
	const uint a0 = (uint)(a), a1 = (uint)(a >> 32);
	const uint b0 = (uint)(b), b1 = (uint)(b >> 32);

	uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c3) : "r" (a1), "r" (b1));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c3) : "r" (c3));

	rax = upsample(c1, c0); rcx = upsample(c3, c2);
#else
        rax = a*b;
        rcx = mul_hi(a,b);
#endif
  
        rax *= Ns;
        rcx += ( (rax != 0)?1:0 );
        rax = mad_hi(rax, N, rcx);

        rcx = rax - N;
        rax = (rax>N)?rcx:rax;

        return rax;
}


/*** Kernel Helpers ***/
// Special thanks to Alex Kruppa for introducing me to Montgomery REDC math!
/* Compute a^{-1} (mod 2^(32 or 64)), according to machine's word size */

ulong invmod2pow_ul (const ulong n)
{
	ulong r;

	const uint in = (uint)n;

	// Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
	// correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
	// Newton iterations are enough.
	r = (n+n+n) ^ ((ulong)2);
	// Newton iteration
	r += r - (ulong)((uint)(r) * (uint)(r) * in);
	r += r - (ulong)((uint)(r) * (uint)(r) * in);
	r += r - (ulong)((uint)(r) * (uint)(r) * in);
	r += r - r * r * n;

	return r;
}


// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
ulong onemod_REDC(const ulong N, ulong rax) {

	ulong rcx;

	// Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
	rcx = (rax!=0)?1:0;
	rax = mad_hi(rax, N, rcx);
	rcx = rax - N;
	rax = (rax>N)?rcx:rax;

	return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
ulong mod_REDC(const ulong a, const ulong N, const ulong Ns) {
	return onemod_REDC(N, Ns*a);
}


// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 6 bits), then just keeps squaring and dividing by 2 when needed.
ulong invpowmod_REDClr (const ulong N, const ulong Ns, const ulong r0, const int bits, const uint nmin) {

	int bbits = bits;
	ulong r = r0;

	// Now work through the other bits of nmin.
	for(; bbits >= 0; --bbits) {
		// Just keep squaring r.
		r = mulmod_REDC(r, r, N, Ns);
		// If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod N).
		if(nmin & (1u << bbits)) {
			r += ( (r&1) ? N : 0 );
			r = r >> 1;
		}
	}

	// Convert back to standard.
	r = mod_REDC (r, N, Ns);

	return r;
}


// Set up to check N's by getting in position with division only.
__kernel void setup(__global ulong * P, __global ulong * Ps, __global ulong * K, __global ulong * lK, __global uint * primecount ) {

	uint gid = get_global_id(0);

	if(gid < primecount[0]){

		ulong my_P = P[gid];

		ulong my_Ps = -invmod2pow_ul (my_P); // Ns = -N^{-1} % 2^64

		// Calculate k0, not in Montgomery form.
		ulong k0 = invpowmod_REDClr(my_P, my_Ps, RS, BBITS, NMIN);

		// calculate k for last value of N, for checksum.
		ulong k1 = invpowmod_REDClr(my_P, my_Ps, RSL, BBITSL, LASTN);

		// store to global arrays
		Ps[gid] = my_Ps;
		K[gid] = k0;
		lK[gid] = k1;

	}

}



