!=======================================================================
MODULE RAY3D_PAR
! This replaces 'ray3d_param.inc'

!
! namelen is the length of filenames (in characters)
! maxlay is the maximum number of layers allowed in the input model
! maxtr is the maximum number of traces allowed
! buffsize is the max. line length assumed for reading files.
      integer namelen, maxlay, maxtr, buffsize
      parameter (namelen=40,maxlay=61,maxtr=200,buffsize=120)
      
! maxlayd is the maximum number of layers used in ray3d routine
      integer maxlayd
      parameter (maxlayd=150)

! units for reading and writing
      integer iounit1,iounit2
      parameter (iounit1=1,iounit2=2)
      
! nsamp is the number of samples per trace.
      integer maxsamp
      !parameter (maxsamp=2500)  !! ray3d
      parameter (maxsamp=2048)  !! Shibutani code

! nsub is the number of sublayers for each horizontal layer
!   to mimic the given velocity gradient
      integer nsub
      parameter (nsub=5)
END MODULE RAY3D_PAR
!=======================================================================
