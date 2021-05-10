FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	f77
FC            =	f77
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	tps

INCS	      =	maxarray.inc mdtop.inc mdvelo.inc pathtopo.inc traject.inc

OBJS	      =	exitt.o             \
		force.o             \
		in_A.o              \
		in_B.o              \
		in_Slice.o          \
		init_path.o         \
		md.o                \
		pathensemble.o      \
		ran_gauss.o         \
		ran_uniform.o       \
		readinput.o         \
		readpath.o          \
		sample_avhb.o       \
		sample_ct.o         \
		sample_dist.o       \
		sample_umbrella.o   \
		shifting.o          \
		shooting.o          \
		sstmm.o             \
		tps.o               \
		writepath.o          

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
