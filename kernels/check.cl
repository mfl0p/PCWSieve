/*

	check kernel

	Bryan Little 2/12/2023

	validates proper operation of the sieve kernel and notifies CPU if there was an error
	also computes the checksum using a local memory reduction

*/


__kernel __attribute__ ((reqd_work_group_size(256, 1, 1))) void check(__global ulong *g_K, __global ulong *g_lK, __global uint *g_primecount, __global ulong *g_P, __global ulong *g_checksum, uint numgroups) {

	uint gid = get_global_id(0);
	uint lid = get_local_id(0);
	__local ulong checksum[256];
	uint pcnt = g_primecount[0];

	if(gid < pcnt){

		ulong my_K = g_K[gid];
		ulong last_K = g_lK[gid];
		ulong my_P = g_P[gid];

		// add my_P and my_K to local memory
		checksum[lid] = my_P + my_K;

		// should match if sieve kernel calculated from nmin to nmax correctly.
		if(my_K != last_K){
			// printf("n %u bbits1 %d r1 %llu checksum mismatch %llu vs %llu\n",my_lastN,bbits1,r1,my_K,kpos);
			// checksum mismatch, set flag
			atomic_or(&g_primecount[3], 1);
		}
	}
	else{
		checksum[lid] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// local memory reduction
	for(int s = 128; s > 0; s >>= 1){
		if(lid < s){
			checksum[lid] += checksum[lid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if(lid == 0){
		uint index = get_group_id(0) + 1;

		if(index < numgroups){	
			// add local checksum to global
			g_checksum[index] += checksum[0];
		}
	}

	if(gid == 0){

		// add primecount to total primecount
		g_checksum[0] += pcnt;

		// store largest kernel prime count
		if( pcnt > g_primecount[1] ){
			g_primecount[1] = pcnt;
		}
	}

}





