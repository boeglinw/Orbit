#------Makefile for orbit205 split into subroutines----------

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

ORBIT205_OBJ = chrien.o thal4.o fcn1.o magfld.o accept.o get_psirel.o Em.o write_orbit.o
INCLUDE = ../include/

F90     = gfortran

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-comments \

#FFLAGS  = -ffixed-line-length-none -fd-lines-as-code \

#FDFLAG = -fd-lines-as-code

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -fno-align-commons -w -fno-automatic  -I$(INCLUDE)

THIS_LIB = liborbit205.a

$(THIS_LIB): $(ORBIT205_OBJ)
	ar rv $(THIS_LIB) $(ORBIT205_OBJ)
	ranlib $(THIS_LIB)

#----------------------------------------------------------
all: $(ORBIT205_OBJ) $(THIS_LIB)

orbit.o: orbit.f

#example of creating object files the long way
#magfld.o: magfld.f
#	$(F90) -c magfld.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(ORBIT205_OBJ)
	rm -f orbit.o
	rm -f $(THIS_LIB)
