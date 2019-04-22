# this file is for general settings such as file list, etc.
# machine-specific settings such as include paths and #defines are in Makefile.local

CC   = mpicc
FC   = mpif90

##CFLAGS +=  -fPIC -Wall -O2 -I$(SRCDIR)
LFLAGS   +=   -fPIC -fopenmp

INC2DFFT += -I/home/qwang/local/2decomp_fft/include
LIB2DFFT += -L/home/qwang/local/2decomp_fft/lib

SRCDIR    = src
INCDIR	  = inc
OBJDIR    = obj
LIBDIR    = lib
EXEDIR	  = ./run/
EXE       = photoNs-lcdm
LIBS	 += -l2decomp_fft  
LIBS	 += -fopenmp  -lm -lpthread

OPTS	 += -O2
OPTS	 += -DLONGSHORT
OPTS	 += -DMYALLTOALLV
OPTS	 += -DPMTHREAD
OPTS	 += -DPERIODIC_CONDITION

#OPTS	 += -Wno-unused-but-set-variable -Wno-unused-variable 
#OPTS	 += -Wno-unused-function -Wno-strict-aliasing -fopenmp

SOURCES	  = photoNs.c domains.c initial.c remotes.c toptree.c \
	    utility.c operator.c snapshot.c partmesh.c  fmm.c
	    
CONV2D	  = conv.f90


OBJECTS	  = $(patsubst %.c, $(OBJDIR)/%.o,$(SOURCES))
OBJFORT	  = $(patsubst %.f90, $(OBJDIR)/%.o,$(CONV2D))

exe: $(OBJECTS) $(OBJFORT)
	@mkdir -p $(EXEDIR)
	$(FC) $(OBJDIR)/*.o $(LIBS) $(LIB2DFFT) -o $(EXEDIR)/$(EXE)

$(OBJDIR)/%.o:  $(SRCDIR)/%.c 
	@mkdir -p $(OBJDIR)
	$(CC) -c $(CFLAGS) $(OPTS) $(INCLUDES) $(LIBS) -I$(INCDIR) -o "$@" "$<"

$(OBJDIR)/%.o:  $(SRCDIR)/%.f90 
	$(FC) -c  $(INC2DFFT) -o "$@" "$<"

.PHONY: demo
demo:
	cd $(EXEDIR); \
	mpirun -np 4 ./$(EXE) ../demo/lcdm_g2.run

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(EXEDIR)/$(EXE)

