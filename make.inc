CC = gcc
CFLAGS = -D_FILE_OFFSET_BITS=64 -pthread -fopenmp -DGOTO -DLINUX

LDFLAGS =
LDLIBS = /home/dfabregat/libs/lapack/lapack_LINUX.a /home/dfabregat/libs/GotoBLAS2/libgoto2.a -lgfortran -gfortranbegin -lpthread -lrt
