#------Makefile for rdpar.f split into subroutines----------

RDPAR_OBJ = rdpar.o rdpar_nml.o
INCLUDE = ../include/

F90     = gfortran

FDFLAG  = -fd-lines-as-comments

#FDFLAG  = -fd-lines-as-code 

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -fno-align-commons -w -fno-automatic -I$(INCLUDE)

THIS_LIB = librdpar.a

$(THIS_LIB): $(RDPAR_OBJ)
	ar rv $(THIS_LIB) $(RDPAR_OBJ)
	ranlib $(THIS_LIB)


#----------------------------------------------------------
all: $(RDPAR_OBJ) $(THIS_LIB)

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

#example of creating object files the long way
#intlim.o: intlim.f
#	$(F90) -c intlim.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(RDPAR_OBJ)
	rm -f $(THIS_LIB)


