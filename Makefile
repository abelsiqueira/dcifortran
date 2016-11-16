# DCI v2.0 8/4/1998

AMD_DIR = $(HOME)/Libraries/AMD

F77 = gfortran
#DCIFLAGS = -g -pg
#DCIFLAGS = -O3 -m64 -parallel -mcmodel=medium -shared-intel -L/home/chico/GotoBLAS -lgoto -liomp5 -lpthread
DCIFLAGS = -O3 -L$(GOTOBLAS) -lpthread -g
BLAS_LIBS = -L$(GOTOBLAS) -lgoto
#BLAS_LIBS = -L/opt/intel/mkl/8.0.1/lib/32 -lmkl_ia32 -lguide -lpthread
AMD_LIB = $(AMD_DIR)/Lib/libamdf77.a
#UMFPACK_LIB = /home/chico/UMFPACKv4.4/UMFPACK/Lib/libumfpack.a

OBJ = dcibfgs.o dciblas.o dcicg.o dcichol.o dcicute.o dcidog.o dci.o dcihoriz.o dciio.o dcisteih.o dcivert.o
# compiling and linking the dcicute program.

all: $(OBJ) dcimain.o
	ar rc libdci.a $(OBJ)
	mkdir -p $(CUTEST)/src/dci
	cp -f dci_cutest.pro $(CUTEST)/packages/$(MYARCH)/double/dci
	chmod a+x $(CUTEST)/packages/$(MYARCH)/double/dci
	cp -f dcimain.f $(CUTEST)/src/dci/dci_main.f
	cp -f dcimain.o $(CUTEST)/objects/$(MYARCH)/double/dci_main.o
	cp -f libdci.a $(CUTEST)/objects/$(MYARCH)/double/

# defining cleaning rule.

#$(AMD_LIB):
	#( cd $(AMD_DIR) ; make fortran )

.PHONY:	clean
clean:
	rm -f *.o *~ core

# defining pattern rules.

%.o:	%.f
	$(F77) $(DCIFLAGS) -c $<
