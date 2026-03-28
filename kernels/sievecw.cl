/*

	sievecw kernel

	Bryan Little Mar 2026
	Ken Brazier August 2010

	sieve for Cullen and Woodall factors

*/

typedef struct {
	ulong p;
	uint n;
	int k;
} factor;

// count trailing zeros
// needed because ctz() is undefined in Nvidia and AMD's CL v1.1 implementation
#define __ctz(_X) \
	31u - clz(_X & -_X)


// 1 if a number mod 15 is not divisible by 2 or 3.
//                           0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
__constant int prime15[] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1 };

bool goodfactor(uint uk, uint n, int c){

	ulong k = uk;

	if(prime15[(uint)(((k<<(n&3))+c)%15)] && (uint)(((k<<(n%3))+c)%7) != 0)
		return true;

	return false;

}

void store_factor(ulong p, uint n, uint k, int c, __global uint *g_primecount,  __global factor *g_factor){
	uint idx = atomic_inc(&g_primecount[2]);
	int sk = (c==1) ? k : -((int)k);
	factor fac = { p, n, sk };
	g_factor[idx] = fac;
}

// For any NSTEP.  not as fast as the 32 and SM versions below due to 64 bit mul

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
// rax is passed in as a * Ns.
ulong shiftmod_REDC (const ulong a, const ulong N, ulong rax)
{
	ulong rcx;

	rax = rax << MONT_NSTEP; // So this is a*Ns*(1<<s) == (a<<s)*Ns.
	rcx = a >> NSTEP;

	rcx += ((rax != 0)?1:0);	// if rax != 0, increase rcx

	rax = mad_hi(rax, N, rcx);

	rcx = rax - N;
	rax = (rax>N)?rcx:rax;

	return rax;
}


__kernel void sievecw(__global ulong *g_P, __global ulong *g_Ps, __global ulong *g_K, __global uint *g_primecount, __global factor *g_factor,
			const uint nstart, const uint nend) {

	uint n = nstart;
	ulong kpos;
	uint i;
	uint gid = get_global_id(0);

	if(gid < g_primecount[0]){
		ulong Ps = g_Ps[gid];
		ulong k0 = g_K[gid];
		ulong my_P = g_P[gid];

		do {
			// Select the even one.
			kpos = (((uint)k0) & 1)?(my_P - k0):k0;

			i = (uint)(kpos);
			if(i != 0){
				i = __ctz(i);
				if(i <= NSTEP){
					if ((((uint)(kpos >> 32))>>i) == 0) {
						uint the_k = (uint)(kpos >> i);
						uint the_n = n + i;
						if(the_k <= the_n){
							while(the_k < the_n){
								the_k <<= 1;
								the_n--;
							}
							if(the_k == the_n && the_n < NMAX) {
								int c = (kpos==k0)?-1:1;
								if( goodfactor(the_k, the_n, c)){
									store_factor(my_P, the_n, the_k, c, g_primecount,  g_factor);
								}
							}
						}
					}
					// if (kpos >> 32))>>i > 0, k is too large.  it cannot be greater than n, which is uint.
				}
			}
			else {
				// if this is called, we already know (uint)(kpos) == 0
				// i is >= 32
				i = (uint)(kpos>>32);
				i = __ctz(i) + 32;
				if(i <= NSTEP){
					uint the_k = (uint)(kpos >> i);
					uint the_n = n + i;
					if(the_k <= the_n){
						while(the_k < the_n){
							the_k <<= 1;
							the_n--;
						}
						if(the_k == the_n && the_n < NMAX) {
							int c = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, c)){
								store_factor(my_P, the_n, the_k, c, g_primecount,  g_factor);
							}
						}
					}
				}
			}

			// Proceed to the K for the next N.
			n += NSTEP;
			k0 = shiftmod_REDC(k0, my_P, k0*Ps);

		} while (n < nend);

		g_K[gid] = k0;  // store k0 to global array
	}
}


// For NSTEP == 32

ulong mad_wide_u32 (const uint a, const uint b, ulong c) {

#ifdef __NV_CL_C_VERSION
	asm volatile ("mad.wide.u32 %0, %1, %2, %0;" : "+l" (c) : "r" (a) , "r" (b));
#else
	c += upsample(mul_hi(a, b), a*b);
#endif

	return c;
}


// Same function, for a constant NSTEP of 32.
ulong shiftmod_REDC32 (ulong rcx, const ulong N, const uint rax)
{
	rcx >>= 32;

	rcx += mad_hi( rax, (uint)N, (uint)((rax!=0)?1:0) );

	rcx = mad_wide_u32((rax),((uint)(N>>32)), rcx);

	rcx = (rcx>N)?(rcx-N):rcx;

	return rcx;
}


__kernel void sievecw32(__global ulong *g_P, __global ulong *g_Ps, __global ulong *g_K, __global uint *g_primecount, __global factor *g_factor,
			const uint nstart, const uint nend) {

	uint i;
	uint n = nstart;
	ulong kpos;
	uint gid = get_global_id(0);

	if(gid < g_primecount[0]){
		ulong Ps = g_Ps[gid];
		ulong k0 = g_K[gid];
		ulong my_P = g_P[gid];
		uint Psh = (uint)Ps;

		do {
			// Select the even one.
			kpos = (((uint)k0) & 1)?(my_P - k0):k0;

			i = (uint)(kpos);
			if(i != 0){
				i = __ctz(i);
				if ((((uint)(kpos >> 32))>>i) == 0) {
					uint the_k = (uint)(kpos >> i);
					uint the_n = n + i;
					if(the_k <= the_n){
						while(the_k < the_n){
							the_k <<= 1;
							the_n--;
						}
						if(the_k == the_n && the_n < NMAX) {
							int c = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, c)){
								store_factor(my_P, the_n, the_k, c, g_primecount,  g_factor);
							}
						}
					}
				}
				// if (kpos >> 32))>>i > 0, k is too large.  it cannot be greater than n, which is uint.
			}
			else {
				// if this is called, we already know (uint)(kpos) == 0
				// i is >= 32, and has to be 32 for this kernel
				uint the_k = (uint)(kpos>>32);
				i = __ctz(the_k) + 32;
				if(i == 32){
					uint the_n = n + 32;
					if(the_k <= the_n){
						while(the_k < the_n){
							the_k <<= 1;
							the_n--;
						}
						if(the_k == the_n && the_n < NMAX) {
							int c = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, c)){
								store_factor(my_P, the_n, the_k, c, g_primecount,  g_factor);
							}
						}
					}
				}
			}

			n += 32;
			k0 = shiftmod_REDC32(k0, my_P, ((uint)k0) * Psh);

		} while (n < nend);

		g_K[gid] = k0;  // store k0 to global array
	}
}


// For NSTEP < 32

// Multiply two 32-bit integers to get a 64-bit result.
ulong mul_wide_u32 (const uint a, const uint b) {

	ulong c;

#ifdef __NV_CL_C_VERSION
	asm volatile ("mul.wide.u32 %0, %1, %2;" : "+l" (c) : "r" (a) , "r" (b));
#else
	c = upsample(mul_hi(a, b), a*b);
#endif

	return c;

}


// Same function for NSTEP < 32. (SMall.)
// Third argument must be passed in as only the low register, as we're effectively left-shifting 32 plus a small number.
ulong shiftmod_REDCsm (ulong rcx, const ulong N, uint rax)
{
	rax <<= SM_MONT_NSTEP;
	rcx >>= NSTEP;
	rcx += (ulong)(mad_hi(rax, (uint)N, (uint)((rax!=0)?1:0) ) );

	rcx += mul_wide_u32(rax, (uint)(N>>32));

	rcx = (rcx>N)?(rcx-N):rcx;
	return rcx;
}


__kernel void sievecwsm(__global ulong *g_P, __global ulong *g_Ps, __global ulong *g_K, __global uint *g_primecount, __global factor *g_factor,
			const uint nstart, const uint nend) {

	uint n = nstart;
	ulong kpos;
	uint i;
	uint gid = get_global_id(0);

	if(gid < g_primecount[0]){
		ulong Ps = g_Ps[gid];
		ulong k0 = g_K[gid];
		ulong my_P = g_P[gid];
		uint Psh = (uint)Ps;

		do {
			// Select the even one.
			kpos = (((uint)k0) & 1)?(my_P - k0):k0;

			i = (uint)(kpos);
			if(i != 0){
				i = __ctz(i);
				if(i <= NSTEP){
					if ((((uint)(kpos >> 32))>>i) == 0) {
						uint the_k = (uint)(kpos >> i);
						uint the_n = n + i;
						if(the_k <= the_n){
							while(the_k < the_n){
								the_k <<= 1;
								the_n--;
							}
							if(the_k == the_n && the_n < NMAX) {
								int c = (kpos==k0)?-1:1;
								if( goodfactor(the_k, the_n, c)){
									store_factor(my_P, the_n, the_k, c, g_primecount,  g_factor);
								}
							}
						}
					}
					// if (kpos >> 32))>>i > 0, k is too large.  it cannot be greater than n, which is uint.
				}
			}
			// if lower 32 bits of kpos are zero, then i will be >= 32 > NSTEP

			n += NSTEP;
			k0 = shiftmod_REDCsm(k0, my_P, ((uint)k0)*Psh);


		} while (n < nend);

		g_K[gid] = k0;  // store k0 to global array
	}
}




