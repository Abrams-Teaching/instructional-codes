FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROG1         =	overdist
PROG2         =	analyse

INCS	      =	maxarray.inc system.inc

OBJS1	      =	disp_mono.o      \
		ener_mono.o      \
		ener_tot.o       \
		grow.o           \
		overdist.o       \
		ran_sphere.o     \
		ran_uniform.o    \
		regrow.o         \
		sstmm.o          \
		widom.o

OBJS2         =	analyse.o

all:		$(PROG1) $(PROG2)

$(PROG1)::	$(OBJS1) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS1) -o $(PROG1)

$(PROG2)::	$(OBJS2) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS2) -o $(PROG2)

clean:;		@rm -f $(OBJS1) $(OBJS2) core
