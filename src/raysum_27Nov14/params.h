!
! namelen is the length of filenames (in characters)
! maxlay is the maximum number of layers allowed in the model
! maxtr is the maximum number of traces allowed
! maxseg: maximum # of segments (should be 3*maxlay for 1st-order
!         multiples
! maxph: maximum number of phases per trace
! buffsize is the max. line length assumed for reading files.
      integer namelen, maxlay, maxtr, maxseg, maxph, buffsize
      parameter (namelen=40,maxlay=60,maxtr=200,maxseg=600)
      parameter (maxph=40000,buffsize=120)
      
! Units for reading and writing
      integer iounit1,iounit2
      parameter (iounit1=1,iounit2=2)
      
! pi: duh. ztol: tolerance considered to be equiv. to zero.
      real pi,ztol
      parameter (pi=3.141592653589793,ztol=1.e-7)
      
! nsamp is the number of samples per trace.
      integer maxsamp
      parameter (maxsamp=2048)
