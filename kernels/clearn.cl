/*

	clearn kernel

	Clears prime counter.

*/


__kernel void clearn(__global uint *g_primecount){

	int i = get_global_id(0);

	if(i==0){
		g_primecount[0]=0; // prime count from prp generator
	}


}



