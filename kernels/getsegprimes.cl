/*

	getsegprimes.cl - Bryan Little 4/2025, montgomery arithmetic by Yves Gallot

	generate a segment of 2-PRPs to test

	Generates a list of base 2 probable primes.  These are "industrial grade primes" requiring ~7 times
 	fewer calculations than testing for primality.  This way we can quickly find candidate primes to
	use for the sieve and remove 2-PRPs later on the CPU.  This compute intensive algorithm is fast on GPU
	when compared to a memory access intensive sieve of Eratosthenes.  Implementing a SoE can require 
	millions of memory accesses that can cause the GPU to stall for hundreds of cycles.

	The following approach is used:

	1) Each thread is two turns of the mod 30 wheel.  This eliminates any numbers divisible by 2, 3, or 5.

	2) Using constant bit sieve arrays the numbers divisible by the small primes from 7 to 113 are removed.
	   A set bit in the array represents a multiple of the prime.  Each thread's mod 30 wheel starting number
	   is modulo the small prime to obtain the correct array index that will tell which of the next 32 odd numbers
	   are divisible by the prime. Since it's a mod 30 wheel only 30 of the positions are used.  The resulting uints
	   are bitwise ORed.  The unset bits in the uint represent numbers that aren't divisible by any of the primes from 7 to 113.

	3) Each thread iterates through it's bitsieve unit using the mod 30 wheel index increment.  If a bit is
	   not set, the number is stored to local memory using a local memory atomic counter.

	4) Packing the numbers in local memory allows all threads to stay busy in the next step, which is performing
	   a base 2 PRP test.  If the number passes the test, it is stored in global memory with an atomic counter along
	   with other constant data that will be used in other kernels.
	
*/

// count trailing zeros long
// needed because ctz() is undefined in Nvidia and AMD's CL v1.1 implementation
#define __ctzl(_X) \
	63u - clz(_X & -_X)

// r0 + 2^64 * r1 = a * b
ulong2 mul_wide(const ulong a, const ulong b){
	ulong2 r;
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
	r.s0 = upsample(c1, c0); r.s1 = upsample(c3, c2);
#else
	r.s0 = a * b; r.s1 = mul_hi(a, b);
#endif
	return r;
}

ulong m_mul(ulong a, ulong b, ulong p, ulong q){
	ulong2 ab = mul_wide(a,b);
	ulong m = ab.s0 * q;
	ulong mp = mul_hi(m,p);
	ulong r = ab.s1 - mp;
	return ( ab.s1 < mp ) ? r + p : r;
}

ulong add(ulong a, ulong b, ulong p){
	ulong r;
	ulong c = (a >= p - b) ? p : 0;
	r = a + b - c;
	return r;
}

ulong invert(ulong p){
	ulong p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

bool strong_prp_two(ulong N){
	ulong q = invert(N);
	ulong one = (-N) % N;
	ulong nmo = N - one;
	ulong two = add(one, one, N);
	int t = __ctzl( (N-1) );
	ulong exp = N >> t;
	ulong curBit = 0x8000000000000000;
	curBit >>= ( clz(exp) + 1 );
	/* If N is prime and N = d*2^t+1, where d is odd, then either
		1.  a^d = 1 (mod N), or
		2.  a^(d*2^s) = -1 (mod N) for some s in 0 <= s < t    */
	ulong a = two;
  	/* r <-- a^d mod N, assuming d odd */
	while( curBit ){
		a = m_mul(a,a,N,q);
		if(exp & curBit){
			a = add(a,a,N);
		}
		curBit >>= 1;
	}
	/* Clause 1. and s = 0 case for clause 2. */
	if(a == one || a == nmo){
		return true;
	}
	/* 0 < s < t cases for clause 2. */
	for(int s = 1; s < t; ++s){
		a = m_mul(a,a,N,q);
		if(a == nmo){
	    		return true;
		}
	}
	return false;
}

// 3 * wheel mod 30
// this way we don't have to check for index wrap around
__constant int wheel[24] = {2, 1, 2, 1, 2, 3, 1, 3, 2, 1, 2, 1, 2, 3, 1, 3, 2, 1, 2, 1, 2, 3, 1, 3};

// bit sieve arrays where set bit represents a multiple of the prime
// arrays represent odd numbers only
__constant uint p7[7] = { 270549121, 2164392968, 135274560, 1082196484, 67637280, 541098242, 33818640 };

__constant uint p11[11] = { 4196353, 134283296, 2098176, 67141648, 2148532736, 33570824, 1074266368, 16785412, 537133184, 8392706, 268566592 };

__constant uint p13[13] = { 67117057, 524352, 33558528, 2147745824, 16779264, 1073872912, 8389632, 536936456, 4194816, 268468228, 2097408, 134234114, 1048704 };

__constant uint p17[17] = { 131073, 33554688, 65536, 16777344, 32768, 8388672, 2147500032, 4194336, 1073750016, 2097168, 536875008, 1048584, 268437504, 524292, 134218752, 262146, 67109376 };

__constant uint p19[19] = { 524289, 268435968, 262144, 134217984, 131072, 67108992, 65536, 33554496, 32768, 16777248, 16384, 8388624, 8192, 4194312, 2147487744, 2097156, 1073743872, 1048578, 536871936 };

__constant uint p23[23] = { 8388609, 2048, 4194304, 1024, 2097152, 512, 1048576, 2147483904, 524288, 1073741952, 262144, 536870976, 131072, 268435488, 65536, 134217744, 32768, 67108872, 16384, 33554436, 8192, 16777218, 4096 };

__constant uint p29[29] = { 536870913, 16384, 268435456, 8192, 134217728, 4096, 67108864, 2048, 33554432, 1024, 16777216, 512, 8388608, 256, 4194304, 128, 2097152, 64, 1048576, 32, 524288, 16, 262144, 8, 131072, 2147483652, 65536, 1073741826, 32768 };

__constant uint p31[31] = { 2147483649, 32768, 1073741824, 16384, 536870912, 8192, 268435456, 4096, 134217728, 2048, 67108864, 1024, 33554432, 512, 16777216, 256, 8388608, 128, 4194304, 64, 2097152, 32, 1048576, 16, 524288, 8, 262144, 4, 131072, 2, 65536 };

__constant uint p37[37] = { 1, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 2147483648, 4096, 1073741824, 2048, 536870912, 1024, 268435456, 512, 134217728, 256, 67108864, 128, 33554432, 64, 16777216, 32, 8388608, 16, 4194304, 8, 2097152, 4, 1048576, 2, 524288 };

__constant uint p41[41] = { 1, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 2147483648, 1024, 1073741824, 512, 536870912, 256, 268435456, 128, 134217728, 64, 67108864, 32, 33554432, 16, 16777216, 8, 8388608, 4, 4194304, 2, 2097152 };

__constant uint p43[43] = { 1, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 2147483648, 512, 1073741824, 256, 536870912, 128, 268435456, 64, 134217728, 32, 67108864, 16, 33554432, 8, 16777216, 4, 8388608, 2, 4194304 };

__constant uint p47[47] = { 1, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 2147483648, 128, 1073741824, 64, 536870912, 32, 268435456, 16, 134217728, 8, 67108864, 4, 33554432, 2, 16777216 };

__constant uint p53[53] = { 1, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 2147483648, 16, 1073741824, 8, 536870912, 4, 268435456, 2, 134217728 };

__constant uint p59[59] = { 1, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 2147483648, 2, 1073741824 };

__constant uint p61[61] = { 1, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 2147483648 };

__constant uint p67[67] = { 1, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p71[71] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p73[73] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p79[79] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p83[83] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p89[89] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p97[97] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p101[101] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p103[103] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p107[107] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p109[109] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__constant uint p113[113] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 1073741824, 0, 536870912, 0, 268435456, 0, 134217728, 0, 67108864, 0, 33554432, 0, 16777216, 0, 8388608, 0, 4194304, 0, 2097152, 0, 1048576, 0, 524288, 0, 262144, 0, 131072, 0, 65536, 0, 32768, 0, 16384, 0, 8192, 0, 4096, 0, 2048, 0, 1024, 0, 512, 0, 256, 0, 128, 0, 64, 0, 32, 0, 16, 0, 8, 0, 4, 0, 2, 0 };

__kernel __attribute__ ((reqd_work_group_size(256, 1, 1))) void getsegprimes(ulong low, ulong high, int wheelidx, __global ulong *g_prime, __global uint *g_primecount){

	const uint gid = get_global_id(0);
	const uint lid = get_local_id(0);
	int idx = wheelidx;
	__local ulong sieved[1900];
	__local int count;

	if(lid == 0){
		count = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// each thread is 2 turns of the mod 30 wheel
	ulong P = low + (gid * 60);
#ifdef CKOVERFLOW
	if(P < low) P = high;
#endif
	ulong end = P + 60;
	if(end > high) end = high;
#ifdef CKOVERFLOW
	if(end < P) end = high;
#endif
	// sieve small primes to 113, this seems optimal
	uint bitsieve = p7[P%7] | p11[P%11] | p13[P%13] | p17[P%17] | p19[P%19] | p23[P%23] | p29[P%29] | p31[P%31]
			| p37[P%37] | p41[P%41] | p43[P%43] | p47[P%47] | p53[P%53] | p59[P%59] | p61[P%61] | p67[P%67]
			| p71[P%71] | p73[P%73] | p79[P%79] | p83[P%83] | p89[P%89] | p97[P%97] | p101[P%101]
			| p103[P%103] | p107[P%107] | p109[P%109] | p113[P%113];

	while(P < end){
		if( (bitsieve & 1) == 0 ){
			sieved[atomic_inc(&count)] = P;
		}

		int inc = wheel[idx++];
		P += inc*2;
		bitsieve >>= inc;
#ifdef CKOVERFLOW
		if(P < low) P = end;
#endif
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for(int pos = lid; pos < count; pos += 256){
		ulong p = sieved[pos];
		if( strong_prp_two(p) ){
			g_prime[ atomic_inc(&g_primecount[0]) ] = p;
		}
	}

}


