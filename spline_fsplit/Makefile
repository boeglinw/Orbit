#------Makefile for spline.f split into subroutines----------

SPLINE_OBJ = banfac.o banslv.o bsplvb.o bspp2d.o eknot.o\
             interv.o ppvalw.o sets2d.o seva2d.o seval.o\
             spl2bc.o spl2pp.o spli2d.o zpline.o spline_rev.o
INCLUDE = ../include/

F90     = gfortran

#FFLAGS = -ffixed-line-length-none -fd-lines-as-comments \

FFLAGS = -ffixed-line-length-none $(FDFLAG) \
         -fno-align-commons -fno-automatic -W -I$(INCLUDE)

THIS_LIB = libspline.a

$(THIS_LIB): $(SPLINE_OBJ)
	ar rv $(THIS_LIB) $(SPLINE_OBJ)
	ranlib $(THIS_LIB)

#----------------------------------------------------------
all: $(SPLINE_OBJ) $(THIS_LIB)

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o  $@

#example of creating object files the long way
#banfac.o: banfac.f
#	$(F90) -c banfac.f $(FFLAGS) -o $@
#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f  $(SPLINE_OBJ)
	rm -f $(THIS_LIB)

