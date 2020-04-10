#------Makefile for plotmed205 split into subroutines----------

PLOTMED205_OBJ = pmedinit.o plotmed2.o plotmed.o zlim.o
INCLUDE = ../include/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -fno-align-commons -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libplotmed205.a

$(THIS_LIB): $(PLOTMED205_OBJ)
	ar rv $(THIS_LIB) $(PLOTMED205_OBJ)
	ranlib $(THIS_LIB)

#----------------------------------------------------------
all: $(PLOTMED205_OBJ) $(THIS_LIB)

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

#example of creating object files the long way
#magfld.o: magfld.f
#	$(F90) -c magfld.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(PLOTMED205_OBJ)
	rm -f $(THIS_LIB)

