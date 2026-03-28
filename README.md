# PCWSieve

PCWSieve by Bryan Little

A BOINC-enabled OpenCL stand-alone sieve for Proth, Riesel, Cullen, and Woodall factors.

+1-1 search algorithm by
* Geoffrey Reynolds
* Ken Brazier, 2009

With contributions by
* Yves Gallot
* Kim Walisch

## Requirements

* OpenCL v1.1
* 64 bit operating system

## How it works

1. Search parameters are given on the command line.
2. A small group of sieve primes are generated on the GPU.
3. The group of primes are tested for factors in the K and N ranges specified.
4. Repeat #2-3 until checkpoint.  Gather factors and checksum data from GPU.
5. Check the factors for validity on the CPU and see if they have any small prime divisiors.
6. Report any factors that pass the CPU tests to factors.txt, along with a checksum at the end.
7. Checksum can be used to compare results in a BOINC quorum.

## Running the program
```
command line options
* -p
* -P		Sieve primes -p <= p < -P < 2^62
* -k
* -K		Sieve for primes k*2^n+/-1 with -k <= k <= -K < 2^31
* -n
* -N		Sieve for primes k*2^n+/-1 with 65 <= -n <= n <= -N < 2^31
* -c		Sieve for Cullen/Woodall factors
* -s or --test	Perform self test to verify proper operation of the program.

Program gets the OpenCL GPU device index from BOINC.  To run stand-alone, the program will
default to GPU 0 unless an init_data.xml is in the same directory with the format:

<app_init_data>
<gpu_type>NVIDIA</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>

or

<app_init_data>
<gpu_type>ATI</gpu_type>
<gpu_device_num>0</gpu_device_num>
</app_init_data>
```

## Related Links

* [PSieve-CUDA](https://github.com/Ken-g6/PSieve-CUDA)
* [primesieve](https://github.com/kimwalisch/primesieve)
