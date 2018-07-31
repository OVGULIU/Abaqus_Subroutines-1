!DEC$ FREEFORM
!============================================================================
! COMMON BLOCK
!============================================================================
!
!----------------------------------------------------------------------------
!
! number of dummy elements for user output
integer, parameter :: NELOUT = 1
! number of variables per element for user output: (9+1+1+1)*KGP
integer, parameter :: NVAROUT = 1 
! number of integration points per element
integer, parameter :: KGP=8
! number of cartesian coordinates per node
integer, parameter :: KCORD=3
! number of nodes per element
integer, parameter :: KNODE=8
! number of slip systems (each system covers both dirrections)
integer, parameter :: KSS = 1
! user varialbe field
double precision :: VARFLD(NVAROUT,NELOUT)
common VARFLD
!
! variation value for numerical tangent (global)
double precision, parameter :: KVNTG = 1.d-08
! abort criterium/tolerance for local iteration
double precision, parameter :: KTOL = 1.d-09
! variation value for numerical tangent (local)
double precision, parameter :: KVNTL = 1.d-07
! local iteration limit
double precision, parameter :: KITER = 200
! shape function matrix (interpolation matrix)
double precision, dimension(KGP,KNODE) :: KNN
common KNN
!
! integration point coordinates
double precision, dimension(KGP,KCORD) :: KGPCORD
common KGPCORD
!
! integration point weights
double precision, dimension(KGP) :: KWT
common KWT
!
! factor for quadrature rule
double precision :: KQUAD
common KQUAD
!
! identity tensor
double precision, dimension(KCORD,KCORD) :: KID
common KID
!
! flag parameter for UVARM routine using MPI
integer :: KFLAG
common KFLAG
!
!----------------------------------------------------------------------------
