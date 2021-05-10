FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	ewald

INCS	      =	system.inc

OBJS	      =	dderfc.o            \
		ewald.o             \
		fourierspace.o      \
		lattice.o           \
		realspace.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
