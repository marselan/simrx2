LINKER=mpiCC
LDFLAGS=-O3 -lgfortran

CXX=mpiCC
CXXFLAGS=-O3

FC=mpiCC
FCFLAGS=-O3

TARGET=simrx

simrx: globalstructs.o PMPenelope2008.o PMPengeom2008.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o configparser.o simrx2.o mainmpi2.o 
	$(LINKER) -o $@ $^ $(LDFLAGS) 

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -f $(TARGET) *.o
