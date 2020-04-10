#------Makefile for bs_ode split into subroutines----------

BSODE_OBJ = bs_ode.o bsstep.o mmid.o pzextr.o rk4.o \
            rkqc.o rzextr.o
INCLUDE = ../include/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -fno-align-commons -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libbsode.a

$(THIS_LIB): $(BSODE_OBJ)
	ar rv $(THIS_LIB) $(BSODE_OBJ)
	ranlib $(THIS_LIB)

#----------------------------------------------------------
all: $(BSODE_OBJ) $(THIS_LIB)

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

#example of creating object files the long way
#magfld.o: magfld.f
#	$(F90) -c magfld.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(BSODE_OBJ)
	rm -f $(THIS_LIB)


