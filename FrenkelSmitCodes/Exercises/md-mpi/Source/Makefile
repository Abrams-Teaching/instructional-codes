FCOMPFLAGS    =	-O -w
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	hf77
FC            =	hf77
CC            = hcc

MAKEFILE      =	Makefile
PROGRAM       =	md_nve

INCS	      =	system.inc

OBJS	      =	init.o              \
		integrate.o         \
		force.o             \
		merge.o             \
		md.o                \
		mdloop.o            \
		ran_uniform.o       \
		rangauss.o          \
		readdat.o           \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM) -lmpi

clean:;		@rm -f $(OBJS) core
