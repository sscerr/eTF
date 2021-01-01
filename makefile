include ./Machine

FILE = eTF.x
# uncomment this for parallel compilation
COMMDIR=comm
# 
# uncomment this for serial compilation
#COMMDIR=comm_serial

MODULES= TF_mod.o

OBJECTS= avanzamento.o BC.o condinit.o cont_n_U.o \
	 corrente.o derx_1.o dery_1.o derx_BC.o \
         griglia.o int_to_char.o \
	 faraday_x.o faraday_y.o faraday_z.o \
	 filtro_open_x.o filtro_per_y.o filtro_per_z.o flr_i.o \
	 init_filtro_open_x.o init_filtro_per.o traspdist.o \
	 inizializza.o moto_xy.o moto_z.o ohm.o outx.o outt.o \
	 parametri.o poisson.o pressure_i.o u_ei.o pressure_e.o \
         trac.o TwoFluids.o pstartup.o 

UTILLIB= Util/ftrout.o Util/tridLU.o
#UTILLIB= Util/ftrout.o

$(FILE) : $(MODULES) $(OBJECTS) COMM UTIL
	$(LNK) $(LNKFLAGS) $(OBJECTS) $(MODULES) -o $(FILE) $(UTILLIB) $(LIBS) $(COMMDIR)/comm.a $(LIB)

COMM:
	( cd $(COMMDIR) ; make )

UTIL:
	( cd Util ; make )


#$(OBJECTS): %.o: %.f90
#	$(FF) -c  $< -O3

veryclean: 
	rm -f $(OBJECTS) $(FILE) *.mod
	( cd $(COMMDIR) ; make clean )
	( cd Util  ; make clean )
	make clean


$(OBJECTS) : $(MODULES)

include ./Rules
