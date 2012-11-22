CFLAGS=-O3 -Wall

# Modify to the location of CVODE library
CPPFLAGS=-I$(HOME)/local/include/
LDFLAGS=-L$(HOME)/local/lib/ -lsundials_cvode -lsundials_nvecserial

concrete-corrosion: concrete-corrosion.c

.PHONY: clean
clean:
	$(RM) concrete-corrosion
