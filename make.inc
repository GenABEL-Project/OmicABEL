CC = gcc
CFLAGS = -static -static-libgcc -D_FILE_OFFSET_BITS=64 -fopenmp -DMKL

LDLIBS = $(MKLROOT)/lib/intel64/librwthmkl.a $(INTELROOT)/compiler/lib/intel64/libiomp5.a
