# DCI v2.0 8/4/1998

AMD_DIR = $(HOME)/Libraries/AMD

F77 = gfortran
#DCIFLAGS = -g -pg
#DCIFLAGS = -O3 -m64 -parallel -mcmodel=medium -shared-intel -L/home/chico/GotoBLAS -lgoto -liomp5 -lpthread
DCIFLAGS = -O3 -L$(GOTOBLAS) -liomp5 -lpthread -g -m32
BLAS_LIBS = -L$(GOTOBLAS) -lgoto
#BLAS_LIBS = -L/opt/intel/mkl/8.0.1/lib/32 -lmkl_ia32 -lguide -lpthread
AMD_LIB = $(AMD_DIR)/Lib/libamdf77.a
#UMFPACK_LIB = /home/chico/UMFPACKv4.4/UMFPACK/Lib/libumfpack.a

OBJ = dcibfgs.o dciblas.o dcicg.o dcichol.o dcicute.o dcidog.o dci.o dcihoriz.o dciio.o dcisteih.o dcivert.o
# compiling and linking the dcicute program.

all: $(OBJ) dcimain.o
	ar rc libdci.a $(OBJ)
	@echo " "
	@echo "CUTEr interaction initiating"
	@echo " "
	mv libdci.a $(MYCUTER)/double/lib
	#cp -f $(AMD_LIB) $(MYCUTER)/double/lib
	mkdir -p $(CUTER)/common/src/pkg/dci
	cp -f dci.sh.pro $(CUTER)/build/prototypes/dci.sh.pro
	sed -f $(MYCUTER)/double/config/script.sed dci.sh.pro > $(MYCUTER)/bin/dci
	chmod a+x $(MYCUTER)/bin/dci
	cp -f dci.spc $(CUTER)/common/src/pkg/dci/dci.spc
	cp -f dcimain.f $(MYCUTER)/double/bin/dcima.f
	cp -f dcimain.o $(MYCUTER)/double/bin/dcima.o
	@echo " "
	@echo "CUTEr interaction complete"
	@echo " "

# defining cleaning rule.

#$(AMD_LIB):
	#( cd $(AMD_DIR) ; make fortran )

.PHONY:	clean
clean:
	rm -f *.o *~ core

# defining pattern rules.

%.o:	%.f
	$(F77) $(DCIFLAGS) -c $<
