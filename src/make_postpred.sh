#!/bin/bash

PASH1=$NVCOMPILERS/Linux_x86_64/22.5/compilers/lib #PATH TO SHARED LIBRARIES
PASH2=/opt/Geopsy.org/lib
PASH3=/opt/Qt/5.14.2/gcc_64/lib

FC="mpifort -C"

$FC -c postpred.f90

$FC -o postpred -Wl,-rpath,$PASH1 -Wl,-rpath,$PASH2 -Wl,-rpath,$PASH3  ./ray3d.o ./ray3d/ray3d_com.o ./ray3d/ray3d_par.o ./ray3d/theo.o ./ray3d/qlayer.o ./obj/alloc_obj.o ./obj/four1.o ./obj/mod_sac_io.o ./obj/nr.o ./obj/nrutil.o ./obj/pythag.o ./obj/read_input.o  ./obj/rjmcmc_com.o ./obj/convlv.o  ./obj/data_type.o  ./obj/fourrow.o ./obj/loglhood.o  ./obj/MT1DforwardB3.o  ./obj/nrtype.o ./obj/quicksort.o  ./obj/realft.o      ./obj/svdcmp.o ./obj/libgpell.a ./obj/libraysum.a  ./obj/libswd.a postpred.o -lstdc++ -LPASH1 -lblas -llapack -L$PASH2 -lQGpCoreTools -lQGpCoreMath -lQGpCoreWave -L$PASH3 -lQt5Gui -lQt5Core -lGL -lpthread 
