#------Makefile for utilities subroutines----------

WLIB_OBJ = wlib.o
INCLUDE = ../include/

F90     = gfortran

FDFLAG  = -fd-lines-as-comments

#FDFLAG  = -fd-lines-as-code 

FFLAGS  = -ffixed-line-length-none $(FDFLAG) \
          -fno-align-commons -w -fno-automatic -I$(INCLUDE)

THIS_LIB = libw.a

$(THIS_LIB): $(WLIB_OBJ)
	ar rv $(THIS_LIB) $(WLIB_OBJ)
	ranlib $(THIS_LIB)


#----------------------------------------------------------
all: $(WLIB_OBJ) $(THIS_LIB)

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

#example of creating object files the long way
#intlim.o: intlim.f
#	$(F90) -c intlim.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f $(WLIB_OBJ)
	rm -f $(THIS_LIB)


