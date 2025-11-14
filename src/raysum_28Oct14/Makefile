## Mac OSX:
#FC = mpif90 -fc=pgfortran -C
FC = mpif90 -fc=pgfortran -fast -O4 -fastsse
#FC = mpif90 -fc=pgfortran -fast -O4 -fastsse -Mfprelaxed -Mipa=fast,inline
## Kispiox:
#FC = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -tp=shanghai-64 -C
#FC = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -tp=shanghai-64 -fast -O4 -fastsse -Mfprelaxed -Mipa=fast,inline

FO=libraysum.a

MODS = 
OBJS = eigenvec.o matrixops.o eispack-cg.o raysum.o readwrite.o \
	phaselist.o buildmodel.o trace.o misfit.o recipes00.o
                  
$(FO): ${MODS} ${OBJS}
	ar urv libraysum.a ${MODS} ${OBJS}
	cp *.a ../obj/

recipes00.o: recipes00.f; $(FC) -c recipes00.f
eigenvec.o: eigenvec.f; $(FC) -c eigenvec.f
matrixops.o: matrixops.f; $(FC) -c matrixops.f
eispack-cg.o: eispack-cg.f; $(FC) -c eispack-cg.f
raysum.o: raysum.f; $(FC) -c raysum.f
readwrite.o: readwrite.f; $(FC) -c readwrite.f
phaselist.o: phaselist.f; $(FC) -c phaselist.f
buildmodel.o: buildmodel.f; $(FC) -c buildmodel.f
trace.o: trace.f; $(FC) -c trace.f
misfit.o: misfit.f; $(FC) -c misfit.f

clean:	
	/bin/rm -f *.o *.mod *.l core *.a
