INTEL_HEADERS=/opt/intel/intel-opencl-1.2-5.0.0.43/opencl-1.2-sdk-5.0.0.43/include
INTEL_LIB=/opt/intel/intel-opencl-1.2-5.0.0.43/opencl/lib64

WINDOWS_HEADERS=C:\Intel\INDE\code_builder_5.1.0.25\include
WINDOWS_LIB=C:\Intel\INDE\code_builder_5.1.0.25\lib\x86

BEIGNET_HEADERS=/usr/local/lib/x86_64-linux-gnu/beignet/include
BEIGNET_LIB=/usr/local/lib/x86_64-linux-gnu/beignet

INTEL_OCL_LIB=OpenCL
BEIGNET_OCL_LIB=cl

# Change here as you prefer
HEADERS=$(BEIGNET_HEADERS)
LIB=$(BEIGNET_LIB)
OCL_LIB=$(BEIGNET_OCL_LIB)

#.PHONY means that those rules are always executed
.PHONY: main debug test all clean


#-g -O0 -s
main:
	gcc $(CFLAGS) $(OPTS) main.cpp fbcl.cpp file_handler.c utils.c -o Polito-SW -I$(HEADERS) -L$(LIB) -l$(OCL_LIB) -lm
	
debug:
	gcc $(CFLAGS) $(OPTS) main.cpp fbcl.cpp file_handler.c utils.c -o Polito-SW -I$(HEADERS) -L$(LIB) -l$(OCL_LIB) -lm -g -O0 -s -D DEBUG=1 -Wall

test:
	gcc $(CFLAGS) test.cpp -o test -I$(HEADERS) -L$(LIB) -l$(OCL_LIB)

file_handler:	
	gcc $(CFLAGS) file_handler.c -o file_handler.out

all:
	gcc $(CFLAGS) main.cpp fbcl.cpp file_handler.c utils.c -o PolitoSW -I$(HEADERS) -Wall -L$(LIB) -l$(OCL_LIB)
clean:
	rm -rf PolitoSW test

