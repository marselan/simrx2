LINKER=mpiCC
LDFLAGS=-O3 -lgfortran

CXX=mpiCC
CXXFLAGS=-O3

FC=mpiCC
FCFLAGS=-O3 -march=native

TARGET=simrx2

simrx: globalstructs.o Penelope.o Pengeom.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o configparser.o simrx2.o mainmpi2.o 
	$(LINKER) -o $@ $^ $(LDFLAGS) 

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(TARGET) *.o
