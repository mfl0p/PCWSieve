CC = g++
LD = $(CC)

.SUFFIXES:
.SUFFIXES: .o .c .h .cl .cpp

VERSION_MAJOR := 2
VERSION_MINOR := 0
date := $(shell powershell.exe get-date -format FileDate)

APP = PCWSieve-win64-v$(VERSION_MAJOR).$(VERSION_MINOR)-$(date).exe

SRC = main.cpp cl_sieve.cpp cl_sieve.h simpleCL.c simpleCL.h kernels/clearn.cl kernels/clearresult.cl kernels/getsegprimes.cl kernels/sieve.cl kernels/sievecw.cl kernels/setup.cl kernels/check.cl factor_proth.c factor_proth.h verify_factor.c verify_factor.h putil.c putil.h mont.c mont.h
KERNEL_HEADERS = kernels/clearn.h kernels/clearresult.h kernels/sieve.h kernels/sievecw.h kernels/setup.h kernels/check.h kernels/getsegprimes.h
OBJ = main.o cl_sieve.o simpleCL.o factor_proth.o verify_factor.o putil.o mont.o

LIBS = OpenCL.dll

BOINC_DIR = C:/mingwbuilds/boinc
BOINC_INC = -I$(BOINC_DIR)/lib -I$(BOINC_DIR)/api -I$(BOINC_DIR) -I$(BOINC_DIR)/win_build
BOINC_LIB = -L$(BOINC_DIR)/lib -L$(BOINC_DIR)/api -L$(BOINC_DIR) -lboinc_opencl -lboinc_api -lboinc

CFLAGS  = -I . -I kernels -O3 -m64 -Wall -DVERSION_MAJOR=\"$(VERSION_MAJOR)\" -DVERSION_MINOR=\"$(VERSION_MINOR)\"
LDFLAGS = $(CFLAGS) -lstdc++ -static

all : clean $(APP)

$(APP) : $(OBJ)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(BOINC_LIB) -o $@ libprimesievewin.a

main.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ main.cpp

cl_sieve.o : $(SRC) $(KERNEL_HEADERS)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ cl_sieve.cpp

factor_proth.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ factor_proth.c

verify_factor.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ verify_factor.c

mont.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ mont.c

simpleCL.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ simpleCL.c

putil.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ putil.c

.cl.h:
	perl cltoh.pl $< > $@

clean :
	del *.o
	del kernels\*.h
	del $(APP)

