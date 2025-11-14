##
## Receiver function inversion makefile 
## Jan Dettmer jand@uvic.ca 2012
##
include ./Makefile.compiler

##############################################################################
##
##  Define binary and object paths:
##
BIN = ./bin
OBJ = ./obj
OBJS = ray3d/ray3d.o ray3d/ray3d_com.o ray3d/ray3d_par.o ray3d/theo.o ray3d/qlayer.o \
       $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/alloc_obj.o $(OBJ)/read_input.o $(OBJ)/loglhood.o \
       $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o $(OBJ)/mod_sac_io.o \
       $(OBJ)/four1.o $(OBJ)/fourrow.o $(OBJ)/realft.o $(OBJ)/convlv.o $(OBJ)/quicksort.o \
       $(OBJ)/pythag.o $(OBJ)/UPDATE_MAPfile.o  $(OBJ)/xcorr.o $(OBJ)/Zxcorr.o \
       $(OBJ)/COVmatEst.o  $(OBJ)/ZCOVmatEst.o $(OBJ)/UpdateCOV.o \
       $(OBJ)/MT1DforwardB3.o $(OBJ)/prjmh_temper_rf.o
OBJSR= ray3d/ray3d.o ray3d/ray3d_com.o ray3d/ray3d_par.o ray3d/theo.o ray3d/qlayer.o \
       $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/alloc_obj.o $(OBJ)/read_input.o $(OBJ)/loglhood.o \
       $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o $(OBJ)/mod_sac_io.o \
       $(OBJ)/four1.o $(OBJ)/fourrow.o $(OBJ)/realft.o $(OBJ)/convlv.o $(OBJ)/quicksort.o \
       $(OBJ)/pythag.o $(OBJ)/replica.o $(OBJ)/UPDATE_MAPfile.o  $(OBJ)/xcorr.o $(OBJ)/Zxcorr.o \
       $(OBJ)/COVmatEst.o  $(OBJ)/ZCOVmatEst.o $(OBJ)/UpdateCOV.o $(OBJ)/MT1DforwardB3.o
FO = $(OBJ)/libraysum.a $(OBJ)/libswd.a
LIB =
##############################################################################

all: bin/prjmh_temper_rf

prjmh_temper_rf: bin/prjmh_temper_rf

replica: bin/replica

clean: clean.done
	rm -f *.done
	(cd ./raysum/;/bin/rm -f *.o *.mod *.l core *.a;)
	(cd ./swd/;/bin/rm -f *.o *.mod *.l core *.a;)
	(cd ./ray3d/ && /bin/rm -f *.o *.mod *.l core *.a;)
clean.done:
	rm -f bin/*
	rm -f *.o *.mod $(OBJ)/*

##############################################################################
#
# Targets for linking of programs:
#
# PRJMH temper Receiver:
#
$(BIN)/prjmh_temper_rf: libraysum libswd $(OBJS) $(FO)
	$(FC)  -o $(BIN)/prjmh_temper_rf $(OBJS) $(FO) $(LIB)

$(BIN)/replica: libraysum libswd $(OBJSR) $(FO)
	$(FC)  -o $(BIN)/replica $(OBJSR) $(FO) $(LIB)

$(OBJ)/prjmh_temper_rf.o: prjmh_temper_rf.f90 
	$(FC) -c -o $(OBJ)/prjmh_temper_rf.o prjmh_temper_rf.f90 $(LIB)

$(OBJ)/replica.o: replica.f90 
	$(FC) -c -o $(OBJ)/replica.o replica.f90 $(LIB)

$(OBJ)/rjmcmc_com.o: rjmcmc_com.f90 $(OBJ)/data_type.o; $(FC) -c rjmcmc_com.f90 -o $(OBJ)/rjmcmc_com.o $(OBJ)/data_type.o $(LIB)

$(OBJ)/data_type.o: data_type.f90; $(FCSERIAL) -c data_type.f90 -o $(OBJ)/data_type.o $(LIB)
	/bin/cp data_type.mod ./ray3d/

ray3d/ray3d_com.o: ray3d/ray3d_com.f90 ray3d/ray3d_par.o $(OBJ)/data_type.o
	(cd ./ray3d/;$(FCSERIAL) -c ray3d_com.f90 -o ray3d_com.o $(OBJ)/data_type.o ray3d_par.o $(LIB);)
	/bin/cp ./ray3d/ray3d_com.mod .

ray3d/ray3d_par.o: ray3d/ray3d_par.f90 $(OBJ)/data_type.o
	(cd ./ray3d/;$(FCSERIAL) -c ray3d_par.f90 -o ray3d_par.o $(OBJ)/data_type.o $(LIB);)
	/bin/cp ./ray3d/ray3d_par.mod .

##############################################################################
#
# Targets for separate subroutines:
#
libraysum:
	(cd ./raysum/;make;)
	/bin/ln -sf ./raysum/*.mod .
libswd:
	(cd ./swd/;make;)
	/bin/ln -sf ./swd/*.mod .

ray3d/ray3d.o: ray3d/ray3d_com.o ray3d/ray3d_par.o
	(cd ./ray3d/;make;)
	/bin/cp ./ray3d/ray3d.o .

ray3d/theo.o: ray3d/ray3d_com.o ray3d/ray3d_par.o
	(cd ./ray3d/;make;)
	/bin/cp ./ray3d/theo.o .

##############################################################################
#
# Targets to compile NR subroutines that are in separate files:
#
$(OBJ)/read_input.o: read_input.f90
	$(FC) -c -o $(OBJ)/read_input.o read_input.f90

$(OBJ)/alloc_obj.o: alloc_obj.f90
	$(FC) -c -o $(OBJ)/alloc_obj.o alloc_obj.f90

$(OBJ)/loglhood.o: $(OBJ)/quicksort.o loglhood.f90
	$(FC) -c -o $(OBJ)/loglhood.o loglhood.f90

$(OBJ)/ludcmp.o: ludcmp.f90
	$(FC) -c -o $(OBJ)/ludcmp.o ludcmp.f90

$(OBJ)/MT1DforwardB3.o: MT1DforwardB3.f90
	$(FC) -c -o $(OBJ)/MT1DforwardB3.o MT1DforwardB3.f90

$(OBJ)/svdcmp.o: svdcmp.f90
	$(FCSERIAL) -c -o $(OBJ)/svdcmp.o svdcmp.f90

$(OBJ)/convlv.o: convlv.f90
	$(FCSERIAL) -c -o $(OBJ)/convlv.o convlv.f90

$(OBJ)/realft.o: realft.f90
	$(FCSERIAL) -c -o $(OBJ)/realft.o realft.f90

$(OBJ)/four1.o: four1.f90
	$(FCSERIAL) -c -o $(OBJ)/four1.o four1.f90

$(OBJ)/fourrow.o: fourrow.f90
	$(FCSERIAL) -c -o $(OBJ)/fourrow.o fourrow.f90

$(OBJ)/pythag.o: pythag.f90
	$(FCSERIAL) -c -o $(OBJ)/pythag.o pythag.f90

$(OBJ)/UPDATE_MAPfile.o: UPDATE_MAPfile.f90
	$(FCSERIAL) -c -o $(OBJ)/UPDATE_MAPfile.o UPDATE_MAPfile.f90

$(OBJ)/xcorr.o: xcorr.f90
	$(FCSERIAL) -c -o $(OBJ)/xcorr.o xcorr.f90

$(OBJ)/Zxcorr.o: Zxcorr.f90
	$(FCSERIAL) -c -o $(OBJ)/Zxcorr.o Zxcorr.f90

$(OBJ)/COVmatEst.o: COVmatEst.f90
	$(FCSERIAL) -c -o $(OBJ)/COVmatEst.o COVmatEst.f90

$(OBJ)/ZCOVmatEst.o: ZCOVmatEst.f90
	$(FCSERIAL) -c -o $(OBJ)/ZCOVmatEst.o ZCOVmatEst.f90

$(OBJ)/UpdateCOV.o: UpdateCOV.f90
	$(FCSERIAL) -c -o $(OBJ)/UpdateCOV.o UpdateCOV.f90

##############################################################################
#
#  Targets for numerical recipes modules:
#
$(OBJ)/nrtype.o: nrtype.f90 
	$(FCSERIAL) -c -o $(OBJ)/nrtype.o nrtype.f90

$(OBJ)/nr.o: nr.f90 $(OBJ)/nrtype.o
	$(FCSERIAL) -c -o $(OBJ)/nr.o nr.f90

$(OBJ)/nrutil.o: nrutil.f90 $(OBJ)/nrtype.o
	$(FCSERIAL) -c -o $(OBJ)/nrutil.o nrutil.f90

$(OBJ)/quicksort.o: quicksort.f90
	$(FCSERIAL) -c -o $(OBJ)/quicksort.o quicksort.f90

$(OBJ)/mod_sac_io.o : mod_sac_io.f90
	$(FCSERIAL) -c -o $(OBJ)/mod_sac_io.o mod_sac_io.f90

##############################################################################
# ...this is the end my fiend
#
# EOF
