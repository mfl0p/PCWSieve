/*

	clearresult kernel

	Bryan Little 2/12/2023

	Clears results

*/


__kernel void clearresult(__global uint *flag, __global uint *factorcount, __global ulong *checksum, __global uint *primecount, uint numgroups){

	int i = get_global_id(0);

	if(i < numgroups){
		checksum[i] = 0;	// index 0 is total primecount between checkpoints.  index 1 to 'numgroups' are for each workgroup's checksum
	}
	
	if(i == 0){
		factorcount[0] = 0;	// # of factors found
		flag[0] = 0;		// set to 1 if there is a gpu checksum error
		primecount[1]=0;	// keep track of largest kernel prime count
	}

}

