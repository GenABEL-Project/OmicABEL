include ./make.inc

SRCDIR = ./src
RESH_SRCDIR = ./src/reshuffle
DRIVER = ./CLAK-GWAS

#QUICK and DIRTY
CXX=g++

CFLAGS+=-g -O2  -Wall -I $(SRCDIR)/  # -D__WORDSIZE=64
CXXFLAGS+=-g -O2  -Wall
LDLIBS += -lm

SRCS = $(SRCDIR)/CLAK_GWAS.c $(SRCDIR)/fgls_chol.c $(SRCDIR)/fgls_eigen.c $(SRCDIR)/wrappers.c $(SRCDIR)/timing.c $(SRCDIR)/statistics.c $(SRCDIR)/REML.c $(SRCDIR)/optimization.c $(SRCDIR)/ooc_BLAS.c $(SRCDIR)/double_buffering.c $(SRCDIR)/utils.c $(SRCDIR)/GWAS.c $(SRCDIR)/databel.c 
OBJS = $(SRCS:.c=.o)
RESH_SRCS=$(RESH_SRCDIR)/main.cpp $(RESH_SRCDIR)/iout_file.cpp $(RESH_SRCDIR)/Parameters.cpp $(RESH_SRCDIR)/reshuffle.cpp $(RESH_SRCDIR)/test.cpp
RESH_OBJS = $(RESH_SRCS:.cpp=.o)

.PHONY: all clean reshuffle

all: $(DRIVER) reshuffle

$(DRIVER): $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

reshuffle: $(RESH_OBJS)
	cd src/reshuffle
	$(CXX) $^ -o $(RESH_SRCDIR)/reshuffle

clean:
	$(RM) $(OBJS)
	$(RM) $(DRIVER) 
	$(RM) $(SRCDIR)/*mod*
	$(RM) $(SRCDIR)/*opari_GPU*
	$(RM) $(RESH_OBJS)
	$(RM) $(RESH_SRCDIR)/reshuffle


src/CLAK_GWAS.o: src/CLAK_GWAS.c src/wrappers.h src/utils.h src/GWAS.h \
 src/databel.h src/timing.h src/REML.h src/fgls_chol.h src/fgls_eigen.h \
 src/double_buffering.h
src/GWAS.o: src/GWAS.c src/utils.h src/GWAS.h src/databel.h src/wrappers.h
src/REML.o: src/REML.c src/options.h src/blas.h src/lapack.h src/ooc_BLAS.h \
 src/wrappers.h src/utils.h src/GWAS.h src/databel.h src/statistics.h \
 src/optimization.h src/REML.h
src/databel.o: src/databel.c src/databel.h src/wrappers.h
src/double_buffering.o: src/double_buffering.c src/GWAS.h src/databel.h \
 src/wrappers.h src/double_buffering.h
src/fgls_chol.o: src/fgls_chol.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/timing.h \
 src/double_buffering.h src/utils.h src/fgls_chol.h
src/fgls_eigen.o: src/fgls_eigen.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/timing.h \
 src/double_buffering.h src/ooc_BLAS.h src/utils.h src/fgls_eigen.h
src/ooc_BLAS.o: src/ooc_BLAS.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/utils.h \
 src/double_buffering.h src/ooc_BLAS.h
src/optimization.o: src/optimization.c src/optimization.h
src/statistics.o: src/statistics.c src/statistics.h
src/timing.o: src/timing.c src/timing.h
src/utils.o: src/utils.c \
 src/GWAS.h src/databel.h src/wrappers.h src/utils.h
src/wrappers.o: src/wrappers.c src/GWAS.h src/databel.h src/wrappers.h
