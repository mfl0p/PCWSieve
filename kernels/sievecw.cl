/*

	sievecw kernel

	Bryan Little 2/12/2023
	Ken Brazier August 2010

	sieve for Cullen and Woodall factors

*/


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


// For any nstep.  not as fast as the 32 and SM versions below

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
// rax is passed in as a * Ns.
ulong shiftmod_REDC (const ulong a, const ulong N, ulong rax, const uint mont_nstep, const uint nstep)
{
	ulong rcx;

	rax = rax << mont_nstep; // So this is a*Ns*(1<<s) == (a<<s)*Ns.
	rcx = a >> nstep;

	rcx += ((rax != 0)?1:0);	// if rax != 0, increase rcx

	rax = mad_hi(rax, N, rcx);

	rcx = rax - N;
	rax = (rax>N)?rcx:rax;

	return rax;
}


__kernel void sievecw(__global ulong * g_P, __global ulong * g_Ps, __global ulong * g_K, __global uint * primecount, __global uint2 * factorKN, __global long * factorP,
			__global uint * factorCnt, const uint N, const uint nstep, const uint kernel_nstep, const uint mont_nstep, const uint nmax, const uint kmin,
			const uint kmax) {

	uint n = N;
	ulong kpos;
	uint i;
	uint l_nmax = n + kernel_nstep;
	if(l_nmax > nmax) l_nmax = nmax;

	uint gid = get_global_id(0);

	if(gid < primecount[0]){
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
				if(i <= nstep){
					if ((((uint)(kpos >> 32))>>i) == 0) {
						uint the_k = (uint)(kpos >> i);
						uint the_n = n + i;
						if(the_k <= the_n){
							while(the_k < the_n){
								the_k <<= 1;
								the_n--;
							}
							if(the_k == the_n && the_n <= l_nmax) {
								int s = (kpos==k0)?-1:1;
								if( goodfactor(the_k, the_n, s)){
									int I = atomic_inc(&factorCnt[0]);
									factorP[I] = (s==1) ? (long)my_P : -((long)my_P);
									factorKN[I] = (uint2){ the_k, the_n };
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
				if(i <= nstep){
					uint the_k = (uint)(kpos >> i);
					uint the_n = n + i;
					if(the_k <= the_n){
						while(the_k < the_n){
							the_k <<= 1;
							the_n--;
						}
						if(the_k == the_n && the_n <= l_nmax) {
							int s = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, s)){
								int I = atomic_inc(&factorCnt[0]);
								factorP[I] = (s==1) ? (long)my_P : -((long)my_P);
								factorKN[I] = (uint2){ the_k, the_n };
							}
						}
					}
				}
			}

			// Proceed to the K for the next N.
			n += nstep;
			k0 = shiftmod_REDC(k0, my_P, ((uint)k0)*Psh, mont_nstep, nstep);

		} while (n < l_nmax);


		g_K[gid] = k0;  // store k0 to global array

	}

}


// For nstep == 32

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


__kernel void sievecw32(__global ulong * g_P, __global ulong * g_Ps, __global ulong * g_K, __global uint * primecount, __global uint2 * factorKN, __global long * factorP,
			__global uint * factorCnt, const uint N, const uint nstep, const uint kernel_nstep, const uint mont_nstep, const uint nmax, const uint kmin,
			const uint kmax) {

	uint i;
	uint n = N;
	ulong kpos;
	uint l_nmax = n + kernel_nstep;
	if(l_nmax > nmax) l_nmax = nmax;

	uint gid = get_global_id(0);

	if(gid < primecount[0]){
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
						if(the_k == the_n && the_n <= l_nmax) {
							int s = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, s)){
								int I = atomic_inc(&factorCnt[0]);
								factorP[I] = (s==1) ? (long)my_P : -((long)my_P);
								factorKN[I] = (uint2){ the_k, the_n };
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
						if(the_k == the_n && the_n <= l_nmax) {
							int s = (kpos==k0)?-1:1;
							if( goodfactor(the_k, the_n, s)){
								int I = atomic_inc(&factorCnt[0]);
								factorP[I] = (s==1) ? (long)my_P : -((long)my_P);
								factorKN[I] = (uint2){ the_k, the_n };
							}
						}
					}
				}
			}

			n += 32;
			k0 = shiftmod_REDC32(k0, my_P, ((uint)k0) * Psh);

		} while (n < l_nmax);


		g_K[gid] = k0;  // store k0 to global array

	}

}


// For nstep < 32




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


// Same function for nstep < 32. (SMall.)
// Third argument must be passed in as only the low register, as we're effectively left-shifting 32 plus a small number.
ulong shiftmod_REDCsm (ulong rcx, const ulong N, uint rax, const uint sm_mont_nstep, const uint nstep)
{
	rax <<= sm_mont_nstep;
	rcx >>= nstep;
	rcx += (ulong)(mad_hi(rax, (uint)N, (uint)((rax!=0)?1:0) ) );

	rcx += mul_wide_u32(rax, (uint)(N>>32));

	rcx = (rcx>N)?(rcx-N):rcx;
	return rcx;
}


__kernel void sievecwsm(__global ulong * g_P, __global ulong * g_Ps, __global ulong * g_K, __global uint * primecount, __global uint2 * factorKN, __global long * factorP,
			__global uint * factorCnt, const uint N, const uint nstep, const uint kernel_nstep, const uint mont_nstep, const uint nmax, const uint kmin,
			const uint kmax) {

	uint n = N;
	ulong kpos;
	uint i;
	const uint sm_mont_nstep = mont_nstep - 32;
	uint l_nmax = n + kernel_nstep;
	if(l_nmax > nmax) l_nmax = nmax;

	uint gid = get_global_id(0);

	if(gid < primecount[0]){
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
				if(i <= nstep){
					if ((((uint)(kpos >> 32))>>i) == 0) {
						uint the_k = (uint)(kpos >> i);
						uint the_n = n + i;
						if(the_k <= the_n){
							while(the_k < the_n){
								the_k <<= 1;
								the_n--;
							}
							if(the_k == the_n && the_n <= l_nmax) {
								int s = (kpos==k0)?-1:1;
								if( goodfactor(the_k, the_n, s)){
									int I = atomic_inc(&factorCnt[0]);
									factorP[I] = (s==1) ? (long)my_P : -((long)my_P);
									factorKN[I] = (uint2){ the_k, the_n };
								}
							}
						}
					}
					// if (kpos >> 32))>>i > 0, k is too large.  it cannot be greater than n, which is uint.
				}
			}
			// if lower 32 bits of kpos are zero, then i will be >= 32 > nstep

			n += nstep;
			k0 = shiftmod_REDCsm(k0, my_P, ((uint)k0)*Psh, sm_mont_nstep, nstep);


		} while (n < l_nmax);


		g_K[gid] = k0;  // store k0 to global array

	}

}




