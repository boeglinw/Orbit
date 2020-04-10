#-----Makefile for orbit205 executable w/o plotting-------

INCLUDE = ./include/

EXEC   = orbit3
EXEC1  = check_field

F90    = gfortran

FDFLAG = -fd-lines-as-comments
export FDFLAG

FFLAGS = -ffixed-line-length-none $(FDFLAG)  \
         -fno-align-commons -fno-automatic -I$(INCLUDE)
LDFLAGS	= -g

SPLINE     = ./spline_fsplit
ORBIT205   = ./orbit205_fsplit
ORBLIM     = ./orblim_fsplit
PLOTMED205 = ./plotmed205_fsplit
RDEQDSK    = ./rdeqdsk_fsplit
RDPAR      = ./rdpar_fsplit
BSODE      = ./bs_ode_fsplit 
UTIL      =  ./utilities

MAIN   =  $(ORBIT205)/orbit.o
MAINSRC   =  $(ORBIT205)/orbit.f

MAIN1   =  $(ORBIT205)/check_field.o
MAINSRC1   =  $(ORBIT205)/check_field.f

FLIBS	= -L$(SPLINE) -L$(RDEQDSK) -L$(ORBLIM) -L$(PLOTMED205) \
          -L$(ORBIT205) -L$(RDPAR) -L$(BSODE) -L$(UTIL)\
          -lorbit205 -lplotmed205 \
          -lorblim -lrdeqdsk -lspline -lbsode -lrdpar -lw

#----------------------------------------------------------

$(EXEC): bsode spline orbit205_f orblim plotmed205 rdeqdsk rdpar util $(MAIN)
	$(F90) $(LDFLAGS) -o ./bin/$(EXEC) $(MAIN) $(FLIBS)


$(EXEC1): bsode spline orbit205_f orblim plotmed205 rdeqdsk rdpar util $(MAIN1)
	$(F90) $(LDFLAGS) -o ./bin/$(EXEC1) $(MAIN1) $(FLIBS)


spline:
	cd spline_fsplit ; make all
orbit205_f:
	cd orbit205_fsplit ; make all
orblim:
	cd orblim_fsplit ; make all
plotmed205:
	cd plotmed205_fsplit ; make all
rdeqdsk:
	cd rdeqdsk_fsplit ; make all
rdpar:
	cd rdpar_fsplit ; make all
bsode: 
	cd bs_ode_fsplit ; make all
util: 
	cd utilities; make all

$(MAIN): $(MAINSRC)
	cd $(ORBIT205); make orbit.o

$(MAIN1): $(MAINSRC1)
	cd $(ORBIT205); make check_field.o

#----------------------------------------------------------
.PHONY : clean

clean: clean_spline clean_orbit205 clean_orblim clean_rdpar \
       clean_plotmed205 clean_rdeqdsk clean_bsode clean_util
	cd bin ; rm -f $(EXEC)

clean_spline:
	cd spline_fsplit ; make clean
clean_orbit205:
	cd orbit205_fsplit ; make clean
clean_orblim:
	cd orblim_fsplit ; make clean
clean_plotmed205:
	cd plotmed205_fsplit ; make clean
clean_rdeqdsk:
	cd rdeqdsk_fsplit ; make clean
clean_rdpar:
	cd rdpar_fsplit ; make clean
clean_bsode: 
	cd bs_ode_fsplit ; make clean
clean_util: 
	cd utilities ; make clean


