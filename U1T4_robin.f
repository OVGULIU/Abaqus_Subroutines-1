!DEC$ FREEFORM
!===============================================================================
! DEVELOPER: EMMA GRIFFITHS
! YEAR: 2017
!===============================================================================

! include of ABAQUS subroutines
INCLUDE 'UEXTERNALDB.f'
!INCLUDE 'UVARM.f'
!INCLUDE 'SDVINI.f'

module functions

    implicit none
    !---------------------------------------------------------------------------          
    public :: sgn,norm,trace,det,inv,dot,ddot,dya,matvec,matmat
    
    contains

    !===========================================================================
    ! sign function of a scalar
    function sgn(sca) 
   
        !-----------------------------------------------------------------------  
        ! declaration

        double precision :: sca  
        double precision :: sgn    
        integer :: i
        !-----------------------------------------------------------------------

        if (sca .gt. 0.d0) then
            sgn = 1.d0
        elseif (sca .lt. 0.d0) then
            sgn = -1.d0
        elseif (sca == 0.d0) then
            sgn = 0.d0
        end if

    end function sgn
    !===========================================================================
    
    !===========================================================================
    ! euclidean norm of a vector
    function norm(vec) 
   
        !-----------------------------------------------------------------------  
        ! declaration

        double precision, dimension(:) :: vec
        double precision :: norm

        integer :: i
        !-----------------------------------------------------------------------

        norm = 0.
        if (size(vec) == 3) then
            norm = (vec(1)**2.d0+vec(2)**2.d0+vec(3)**2.d0)**(0.5d0)
        elseif (size(vec) == 2) then
            norm = (vec(1)**2.d0+vec(2)**2.d0)**(0.5d0)
        else
            do i=1,size(vec)
                norm = norm + vec(i)**2.d0
            end do
            norm = sqrt(norm)
        end if

    end function norm
   !============================================================================

    !===========================================================================
    ! trace of a tensor
    function trace(mat) 
   
        !-----------------------------------------------------------------------  
        ! declaration

        double precision, dimension(:,:) :: mat  
        double precision :: trace    
        integer :: i
        !-----------------------------------------------------------------------

        trace = 0.d0
        if (size(mat,1) == size(mat,2)) then
            do i=1,size(mat,1)
                trace = trace + mat(i,i)
            end do
        else
            stop "Error in function `trace' - matrix is non-sqare!"
        end if

    end function trace
    !===========================================================================

    !===========================================================================
    ! determinant of a tensor
    function det(mat)

        !-----------------------------------------------------------------------    
        ! declaration

        double precision, dimension(:,:) :: mat
        double precision :: det
        !-----------------------------------------------------------------------

        det = 0.d0
        if (size(mat,1)==size(mat,2) .and. size(mat,1)==3) then
            det = mat(1,1)*mat(2,2)*mat(3,3) &
                + mat(1,2)*mat(2,3)*mat(3,1) &
                + mat(1,3)*mat(2,1)*mat(3,2) &
                - mat(1,3)*mat(2,2)*mat(3,1) &
                - mat(1,2)*mat(2,1)*mat(3,3) &
                - mat(1,1)*mat(2,3)*mat(3,2)
        else
            stop "Error in function `det' - tensor is non-sqare!"
        end if
            
    end function det
    !===========================================================================

    !===========================================================================
    ! inverse of a tensor
    function inv(mat)
   
        !-----------------------------------------------------------------------    
        ! declaration

        double precision, dimension(:,:) :: mat
        double precision, dimension(size(mat,1),size(mat,2)) :: inv
        !-----------------------------------------------------------------------

        inv = 0.d0
        if (size(mat,1)==size(mat,2) .and. size(mat,1)==3 ) then

            inv(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
            inv(1,2) = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
            inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
        
            inv(2,1) = mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
            inv(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
            inv(2,3) = mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
        
            inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
            inv(3,2) = mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)
            inv(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
        
            inv=inv/det(mat)
            
        elseif (size(mat,1)==size(mat,2) .and. size(mat,1)==2 ) then
            
            inv(1,1) =  mat(2,2)
            inv(1,2) = -mat(1,2)
            inv(2,1) = -mat(2,1)
            inv(2,2) =  mat(1,1)
            
            inv=inv/det(mat)
            
        elseif (size(mat,1)==size(mat,2) .and. size(mat,1)==1 ) then
            
            inv(1,1) = 1.d0/mat(1,1)
            
        else
            stop "Error in function `inv' - tensor is non-sqare or larger then 3x3!"
        end if

    end function inv
    !===========================================================================

    !===========================================================================
    ! scalar/dot product of two vectors of the same dimension
    function dot(vec1,vec2)

        !-----------------------------------------------------------------------   
        ! declaration

        double precision, dimension(:) :: vec1,vec2
        double precision :: dot

        integer :: i
        !-----------------------------------------------------------------------

        dot = 0.d0
        if (size(vec1)==size(vec2)) then  
            do i=1,size(vec1)
                dot = dot + vec1(i)*vec2(i)
            end do
        else
            stop "Error in function `dot' - vectors have not the same length!"
        end if

    end function dot
    !===========================================================================

    !===========================================================================
    ! scalar/dot product of two tensors of the same dimensions
    function ddot(mat1,mat2)

        !-----------------------------------------------------------------------   
        ! declaration

        double precision, dimension(:,:) :: mat1,mat2
        double precision :: ddot

        integer :: i,j
        !-----------------------------------------------------------------------

        ddot = 0.d0
        if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2)) then  
            do i=1,size(mat1,1)
                do j=1,size(mat1,2)
                    ddot = ddot + mat1(i,j)*mat2(i,j)
                end do
            end do
        else
            stop "Error in function `ddot' - tensor dimensions are not the same!"
        end if

    end function ddot
    !===========================================================================

    !===========================================================================
    ! dyadic prodcut of two vectors of the same length
    function dya(vec1,vec2)
   
        !-----------------------------------------------------------------------    
        ! declaration  

        integer :: i,j   
        double precision, dimension(:) :: vec1,vec2
        double precision, dimension(size(vec1),size(vec2)) :: dya
        !-----------------------------------------------------------------------

        dya = 0.d0
        if (size(vec1)==size(vec2)) then

            do i=1,size(vec1)
                do j=1,size(vec1)
    
                    dya(i,j) = vec1(i)*vec2(j)
    
                end do
            end do

        else
            stop "Error in function `dya' - vector lengths are not the same!"
        end if

    end function dya
    !===========================================================================
    
    !============================================================================
    ! cross product of two vectors of the same dimension
    function cross(vec1,vec2)

        !-----------------------------------------------------------------------   
        ! declaration

        double precision, dimension(:) :: vec1,vec2
        double precision, dimension(size(vec1)) :: cross

        integer :: i
        !-----------------------------------------------------------------------

        cross = 0.d0
        if ((size(vec1)==size(vec2)) .and. (size(vec1)==3)) then  
            cross(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
            cross(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
            cross(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
        else
            stop "error in `cross' - vector lengths are not the same or not given in 3D!"
        end if

    end function cross
    !===========================================================================

    !===========================================================================
    ! matrix-vector operation
    function matvec(mat,vec)
   
        !-----------------------------------------------------------------------    
        ! declaration

        double precision, dimension(:,:) :: mat
        double precision, dimension(:) :: vec
        double precision, allocatable, dimension(:) :: matvec      
        integer :: i,j
        !-----------------------------------------------------------------------

        if (size(mat,2) == size(vec)) then

            allocate(matvec(size(mat,1)))
            matvec = 0.d0
            do i=1,size(mat,1)
                do j=1,size(mat,2)
                    matvec(i) = matvec(i) + mat(i,j)*vec(j)
                end do
            end do
    
        else 
            stop "Dimension error in function `matvec' - dimension of vector must be consistent with the dimensions of the matrix"
        end if
 
    end function matvec
    !===========================================================================

    !===========================================================================
    ! matrix-matrix operation (only of square matrices)
    function matmat(mat1,mat2)
   
        !-----------------------------------------------------------------------    
        ! declaration

        double precision, dimension(:,:) :: mat1
        double precision, dimension(:,:) :: mat2        
        double precision, dimension(size(mat1,1),size(mat2,2)) :: matmat
        !-----------------------------------------------------------------------

        matmat=0.d0
        if (size(mat1,2) == size(mat2,1)) then

            matmat(1,1) = mat1(1,1)*mat2(1,1)+mat1(1,2)*mat2(2,1)+mat1(1,3)*mat2(3,1)
            matmat(1,2) = mat1(1,1)*mat2(1,2)+mat1(1,2)*mat2(2,2)+mat1(1,3)*mat2(3,2)
            matmat(1,3) = mat1(1,1)*mat2(1,3)+mat1(1,2)*mat2(2,3)+mat1(1,3)*mat2(3,3)
        
            matmat(2,1) = mat1(2,1)*mat2(1,1)+mat1(2,2)*mat2(2,1)+mat1(2,3)*mat2(3,1)
            matmat(2,2) = mat1(2,1)*mat2(1,2)+mat1(2,2)*mat2(2,2)+mat1(2,3)*mat2(3,2)
            matmat(2,3) = mat1(2,1)*mat2(1,3)+mat1(2,2)*mat2(2,3)+mat1(2,3)*mat2(3,3)
        
            matmat(3,1) = mat1(3,1)*mat2(1,1)+mat1(3,2)*mat2(2,1)+mat1(3,3)*mat2(3,1)
            matmat(3,2) = mat1(3,1)*mat2(1,2)+mat1(3,2)*mat2(2,2)+mat1(3,3)*mat2(3,2)
            matmat(3,3) = mat1(3,1)*mat2(1,3)+mat1(3,2)*mat2(2,3)+mat1(3,3)*mat2(3,3)

        else 
            stop "Dimension error in function `matmat' - matrix dimensions are not consistent or larger than 3x3"
        end if

    end function matmat
    !===========================================================================

end module functions

!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
! USER SUBROUTINE - USER ELEMENT

! import problem values from abaqus
SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
               PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME, &
               KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF, &
               LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

    use functions

    implicit none
    !-------------------------------------------------------------------------------   
    ! declaration

    ! Abaqus variables
    ! ================  
    ! parameters - numbers
    double precision, parameter :: zero=0.d0
    double precision, parameter :: half=0.5d0
    double precision, parameter :: one=1.d0
    double precision, parameter :: two=2.d0
    double precision, parameter :: three=3.d0
    double precision, parameter :: four=4.d0
    double precision, parameter :: six=6.d0

    ! parameters - problem specification
    integer, parameter :: iCORD=3	!Degrees of freedom (mechanical)
    integer, parameter :: ICORDTOTAL=4  !Degrees of freedom (total per node (3 mech; 1 temp))
    integer, parameter :: iNODE=4	!Number of nodes per element

    ! variables passed in for information (scalar parameters)
    integer, intent(in ) :: NDOFEL         ! number of dof in the element
    integer, intent(in ) :: MLVARX         ! dimensioning parameter ???
    double precision, intent(in ) :: DTIME ! current time increment
    integer, intent(in ) :: PERIOD 	       ! time period of current step
    integer, intent(in ) :: NRHS           ! number of load vectors
    integer, intent(in ) :: NSVARS 	       ! user-defined number of solution-dependent 
                        ! state variables
    integer, intent(in ) :: NPROPS      ! user-defined number of real property 
                        ! values
    integer, intent(in ) :: NJPROP      ! user-defined number of integer property 
                                        ! values 
    integer, intent(in ) :: MCRD        ! max. number of user-defined coordinates at 
                        ! each node point
    integer, intent(in ) :: NNODE       ! user-defined number of nodes on the 
                        ! element
    integer, intent(in ) :: JTYPE       ! integer defining the element type
    integer, intent(in ) :: KSTEP 	    ! current step number:
    integer, intent(in ) :: KINC 	    ! current increment number
    integer, intent(in ) :: JELEM 	    ! element number
    integer, intent(in ) :: NDLOAD 	    ! identification number of load active on 
                        ! this element
    integer, intent(in ) :: MDLOAD      ! total number of loads defined on this 
                        ! element
    integer, intent(in ) :: NPREDF      ! number of predefined field variables

    ! variables to be defined
    double precision, dimension(MLVARX) :: RHS 	     ! residual force vector 
    double precision, dimension(MLVARX) :: RHS_test 	     ! residual force vector 
    double precision, dimension(NDOFEL,NDOFEL) :: AMATRX   ! element stiffness matrix 
    double precision, dimension(NSVARS) :: SVARS  	! values of solution-depending 
                                                    ! state variables
    double precision, dimension(8) :: ENERGY		! element energy quantities 

    ! variables which can be updated
    double precision :: PNEWDT	    ! time step control variable

    ! variables passed in for information (arrays)
    double precision, dimension(NPROPS) :: PROPS      ! real property values
    integer, dimension(NJPROP) :: JPROPS            ! integer property values
    double precision, dimension(MCRD,NNODE) :: COORDS ! original coordinates 
                                                    ! (reference configuration)
    double precision, dimension(NDOFEL) :: U  ! solution variables (displacement, 
                                            ! temperature, ...)
    double precision, dimension(MLVARX) :: DU ! increment values
    double precision, dimension(NDOFEL) :: V  ! time rate change of solution variables
    double precision, dimension(NDOFEL) :: A  ! accelerations of the variables
    double precision, dimension(MDLOAD) :: JDLTYP ! load values
    double precision, dimension(MDLOAD) :: ADLMAG ! total load magnitude
    double precision, dimension(MDLOAD) :: DDLMAG ! increments in the load magnitud
    double precision, dimension(2,NPREDF,NNODE) :: PREDEF ! array of predefined field
                                                    ! variables (nodal temperature)
    integer, dimension(3) :: PARAMS  ! parameter of solution precedure
    integer, dimension(5) :: LFLAGS  ! flags for current solution precedure
                                            ! and requirements for element 
                                            ! calculation
    double precision, dimension(2) :: TIME    ! TIME(1): current value of step time
                                            ! TIME(2): current value of total time

    ! internal variables
    ! ==================

    ! include others
    INCLUDE 'COMMON.f'

    !===============================================================================

    ! FE variables
    double precision :: xi1,xi2,xi3     ! natural coordinates
    double precision :: X1(KNODE)       ! physical coordinates
    double precision :: X2(KNODE)       ! physical coordinates
    double precision :: X3(KNODE)       ! physical coordinates
    double precision :: dNdXi1(KNODE)   ! shape function derivatives
    double precision :: dNdXi2(KNODE)   ! shape function derivatives
    double precision :: dNdXi3(KNODE)   ! shape function derivatives
    double precision :: dNdX1(KGP,KNODE)! shape function derivatives
    double precision :: dNdX2(KGP,KNODE)! shape function derivatives
    double precision :: dNdX3(KGP,KNODE)! shape function derivatives
    double precision :: dX1dxi1         ! derivatives
    double precision :: dX1dxi2         ! derivatives
    double precision :: dX1dxi3         ! derivatives
    double precision :: dX2dxi1         ! derivatives
    double precision :: dX2dxi2         ! derivatives
    double precision :: dX2dxi3         ! derivatives
    double precision :: dX3dxi1         ! derivatives
    double precision :: dX3dxi2         ! derivatives
    double precision :: dX3dxi3         ! derivatives
    double precision :: JJ(KCORD,KCORD) ! Jacobi matrix
    double precision :: detJ(KGP)       ! Jacobi-determinant (reference)
    double precision :: kPhi(KNODE)     ! Electric Potential
    double precision :: CoNODE(KNODE)	! Concentration (from predefined field)
    double precision :: E(KCORD)      	! Electric Field
    double precision :: D(KCORD)      	! Electric Displacement
    double precision :: Qf     	! Free Charge Density
    double precision :: vMeq
    double precision :: dofni(KNODE)        ! current dof
    double precision :: kCo ! Concentration
    double precision :: kEpZero, kEpR, kF, kZ, ksat ! Properties included from input file

    ! trash variables for numerical tangent calculation
    double precision, dimension(MLVARX) :: RHS_num   ! variated residual vector 
    double precision, dimension(NDOFEL) :: U_num     ! variated unknown vector 
    double precision, dimension(NDOFEL) :: DU_num    ! variated increment vector


    ! integer
    integer :: ip,nn,ni,i,dof,K1,K2,j

    integer :: CordRange(iCORD) = (/(i, i=1,3, 1)/)
    integer :: NODERange(iNODE) = (/(i, i=1,4, 1)/)
    
    ! parameters
    double precision, parameter :: VPGP = KCORD*KCORD

    ! ==================
    ! Variables for Robin Boundary Conditions
    ! ==================

    ! Variables for normal calculations
    double precision :: AB_edge(iCORD)           ! Edge between nodes 1 and 2
    double precision :: BC_edge(iCORD)           ! Edge between nodes 2 and 3
    double precision :: CD_edge(iCORD)           ! Edge between nodes 3 and 4
    double precision :: AC_edge(iCORD)           ! Edge between nodes 1 and 3
    double precision :: BD_edge(iCORD)           ! Edge between nodes 2 and 4

    double precision :: N1(iCORD)           ! Normal described by face given by vertex: 1,2,3
    double precision :: N2(iCORD)           ! Normal described by face given by vertex: 2,3,4
    double precision :: N3(iCORD)           ! Normal described by face given by vertex: 1,3,4
    double precision :: N4(iCORD)           ! Normal described by face given by vertex: 1,2,4
    double precision :: N_all(iCORD)        ! Normals described by different faces
    
    double precision :: NODES(iCORD)           ! Nodes of face upon which traction is applied

    integer, parameter :: iGPtri=3
    ! user variables- For triangular elements
    ! ==================
    
    double precision :: pGPCORDtri(iGPtri,iCORD-1)  	! integration point coordinates
    double precision :: pQUADtri               		! factor for quadrature rule
    double precision :: pWTtri(iGPtri)            	! integration point weights
    double precision :: pNNtri(iGPtri,iNODE-1)      	! shape function matrix (interpolation matrix)

    double precision :: xi1tri,xi2tri,xi3tri         ! natural coordinates
    double precision :: DetjTri           ! Jacobi-determinant (reference)
    double precision :: pRobin		! Robin boundary parameter
    ! integer
    integer :: nj,nk,nl,iptri, Filesize 

    ! Allocatable arrays
    integer, DIMENSION(:), ALLOCATABLE :: FrontEle 
    integer, DIMENSION(:), ALLOCATABLE :: BackEle

    !-------------------------------------------------------------------------------  

    ! Properties
    kEpZero = props(1)
    kEpR = props(2)
    kF = props(3)
    
    kZ = props(4)
    kSat = props(5)
    
    pRobin = props(6)
    
    ! element ingoing note =======================
    if (kinc == 1) then
!          	write(*,*)"ndofel ",NDOFEL
    end if
        ! ============================================
        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetFront.inp',status='old')!
    READ(107,*) Filesize
    Allocate ( FrontEle(Filesize) )
    close(107)
    
    open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetFront.csv',status='old')!
    READ(107,*) FrontEle
    close(107)
    
    open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetBack.inp',status='old')!
    READ(107,*) Filesize
    Allocate ( BackEle(Filesize) )
    close(107)
    
    open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetBack.csv',status='old')!
    READ(107,*) BackEle
    close(107)
    ! loop over all integration points (computation of FE variables)
    do ip=1,KGP ! ----------------------------------------------------------
        
        ! get solution-dependent state variables (history variables)
        
        ! natural coordinates of current ip
        xi1 = KGPCORD(ip,1)
        xi2 = KGPCORD(ip,2)
        xi3 = KGPCORD(ip,3)
        
        ! coordinate vectors of current element
        X1 = COORDS(1,:)
        X2 = COORDS(2,:)
        x3 = COORDS(3,:)
    
        ! ------------------------------------------------------------------
        if (KNODE==4) then
    
            ! derivatives of shape functions with respect to natural coordinates                
            dNdXi1(1) = -1.d0
            dNdXi2(1) = -1.d0
            dNdXi3(1) = -1.d0
        
            dNdXi1(2) =  1.d0
            dNdXi2(2) =  0.d0
            dNdXi3(2) =  0.d0
        
            dNdXi1(3) =  0.d0
            dNdXi2(3) =  1.d0
            dNdXi3(3) =  0.d0
        
            dNdXi1(4) =  0.d0
            dNdXi2(4) =  0.d0
            dNdXi3(4) =  1.d0
                     
        else 
            write(*,*) "Error in computation of shape function derivatives. The number of nodes does not conform with the element type (4 node tetrahedral element)."  
            call XIT   
        end if
    
        ! derivatives of physical coordinates with respect to natural coordinates                
        dX1dxi1=dot(X1,dNdXi1)
        dX1dxi2=dot(X1,dNdXi2)
        dX1dxi3=dot(X1,dNdXi3)
    
        dX2dxi1=dot(X2,dNdXi1)
        dX2dxi2=dot(X2,dNdXi2)
        dX2dxi3=dot(X2,dNdXi3)
    
        dX3dxi1=dot(X3,dNdXi1)
        dX3dxi2=dot(X3,dNdXi2)
        dX3dxi3=dot(X3,dNdXi3)
    
        ! Jacobian determinant (detJ = 6V)
        detJ(ip) = dX1dxi1*dX2dxi2*dX3dxi3 + dX2dxi1*dX3dxi2*dX1dxi3 + dX3dxi1*dX1dxi2*dX2dxi3 &
                 - dX1dxi3*dX2dxi2*dX3dxi1 - dX2dxi3*dX3dxi2*dX1dxi1 - dX3dxi3*dX1dxi2*dX2dxi1
    
        ! derivatives of shape functions with respect to physical coordinates  
        do nn=1,KNODE
            dNdX1(ip,nn) = 1.d0/detJ(ip)*( (dX2dxi2*dX3dxi3-dX3dxi2*dX2dxi3)*dNdXi1(nn) &
                                         + (dX3dxi1*dX2dxi3-dX2dxi1*dX3dxi3)*dNdXi2(nn) &
                                         + (dX2dxi1*dX3dxi2-dX3dxi1*dX2dxi2)*dNdXi3(nn) )
            dNdX2(ip,nn) = 1.d0/detJ(ip)*( (dX3dxi2*dX1dxi3-dX1dxi2*dX3dxi3)*dNdXi1(nn) &
                                         + (dX1dxi1*dX3dxi3-dX3dxi1*dX1dxi3)*dNdXi2(nn) &
                                         + (dX3dxi1*dX1dxi2-dX1dxi1*dX3dxi2)*dNdXi3(nn) )
            dNdX3(ip,nn) = 1.d0/detJ(ip)*( (dX1dxi2*dX2dxi3-dX2dxi2*dX1dxi3)*dNdXi1(nn) &
                                         + (dX2dxi1*dX1dxi3-dX1dxi1*dX2dxi3)*dNdXi2(nn) &
                                         + (dX1dxi1*dX2dxi2-dX2dxi1*dX1dxi2)*dNdXi3(nn) )
        end do

    end do ! ---------------------------------------------------------------

    ! loop over all degrees of freedom
    !if (Lflags(3).eq.5 .OR. Lflags(3).eq.1) then
    RHS=0.d0
    do ip=1,KGP ! ----------------------------------------------------------

        ! Electric Field
        E = 0.d0
        do ni=1,KNODE
            E = E - U(ni)*( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
        end do
!		WRITE(*,*) E
        CoNODE = predef(1,2,:)
        kCo = dot(kNN(ip,:),CoNODE)
        Qf = kF*((kZ*kCo)+(kSat*1.d0))
!		Qf = kF*((kZ*CoNODE)+( kSat*(/1.d0,1.d0,1.d0,1.d0/) ))

        ! Electric displacement field
        D = kEpZero*kEpR*E 
!		write(*,*) "D", D
    ! Free charge density
        if (jElem==41151 .OR. jElem==46798) then
!			open(unit=103, file='/home/cerecam/Desktop/GPG_Cube/Resultscheck.txt', status="old", position="append", action="write")!	
!			write(*,*) "jElem", jElem
!			write(103,*) "jElem", jElem
!			write(103,*) "U", U	
!			write(103,*) "predef(1,2,:)", predef(1,2,:)	
!!			write(103,*) "Qf", Qf
!			write(103,*) "D", D
!			write(103,*) ""
!			close(103)
        end if		

        do ni=1,KNODE !-----------------------------loop-i------------------
            ! internal residual force vector
!           RHS(ni) = RHS(ni) - (KQUAD*KWT(ip)*detJ(ip)*( Qf*kNN(ip,ni)))
            RHS(ni) = RHS(ni) + KQUAD*KWT(ip)*detJ(ip)*dot(D,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))+KQUAD*KWT(ip)*detJ(ip)*Qf*kNN(ip,ni)
!           RHS(ni) = RHS(ni) + KQUAD*KWT(ip)*detJ(ip)*dot(D,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))
!			if (jelem==1) then
!				write(*,*) "ERROR: RHS", dot(D,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))
!			end if

        end do !------------------------------end-loop-ni--------------------
!		if (jElem==56618) then	
!			write(*,*) "RHS", RHS
!		end if
    end do ! ----------------------------end-loop-ip------------------------
            
        
    !if (Lflags(3).eq.2 .OR. Lflags(3).eq.1) then
    RHS_test=0.d0
    do ip=1,KGP ! ----------------------------------------------------------
                
    ! Electric Field
        E = 0.d0
        do ni=1,KNODE
            E = E - U(ni)*( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
        end do

        D = kEpZero*kEpR*E

        do ni=1,KNODE !-----------------------------loop-i------------------
            ! internal residual force vector
            RHS_test(ni) = RHS_test(ni) + KQUAD*KWT(ip)*detJ(ip)*dot(D,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))

        end do !------------------------------end-loop-i--------------------
!		if (jElem==56618) then	
!			write(*,*) "RHS_test", RHS_test	
!		end if
    end do ! ----------------------------end-loop-ip------------------------
! 	AMATRX=0.d0
    do dof=1,NDOFEL ! ------------------------------------------------------
    ! loop over all integration points (computation of residuum)                
        RHS_num=0.d0
        U_NUM = U(1:NDOFEL)
        U_NUM(dof) = U_NUM(dof)+KVNTG
!                DU_NUM = U_NUM(U(1:NDOFEL)-DU(1:NDOFEL))
        ! loop over all integration points (computation of disturbed residuum)
        do ip=1,KGP ! ------------------------------------------------------

            ! Electric Field
            E = 0.d0
            do ni=1,KNODE
                E = E - U_NUM(ni)*( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
            end do
            ! Electric displacement field
            D = kEpZero*kEpR*E

            ! summation over node_i
            do ni=1,KNODE !-----------------------------loop-i------------------
                ! internal residual force vector
                RHS_num(ni) = RHS_num(ni) + KQUAD*KWT(ip)*detJ(ip)*dot(D,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))

            end do !------------------------------end-loop-i--------------------
    
        end do ! ----------------------------end-loop-ip------------------------

        if ((RHS_num(1)-RHS_test(1)).ne.0) then
        !			write(*,*) "RHS_test:	",RHS_test(1:NDOFEL)
        !			write(*,*) "RHS_num:	",RHS_num(1:NDOFEL)
        !			write(*,*) "RHS_diff:	",RHS_num(1:NDOFEL)-RHS_test(1:NDOFEL)
        end if

       AMATRX(:,dof) = - (RHS_num(1:NDOFEL)-RHS_test(1:NDOFEL))/KVNTG

    end do ! ----------------------------end-loop-dof---------------------------
    if (ANY(jElem.eq.BackEle)) then
!		if (ANY(jElem.eq.FrontEle)) then
!			if (ANY(jElem.eq.(/26528,26271,26542,26421,26327,26338,26396,26167,26301,26111,26311,26493,26454,26413,26492,26413,26489,26485,26372,26088,26379,26270,26403,26261,26121,26343, 26148,26333,26516/))) then
!			if (jElem.eq.26516) then
        pGPCORDtri(1,:) = (/ one/six, one/six /)
        pGPCORDtri(2,:) = (/ two/three, one/six /)
        pGPCORDtri(3,:) = (/ one/six, two/three /)
                
                    
        pWTtri = (/ one/six, one/six, one/six/)
    
        do iptri=1,iGPtri
            xi1tri=pGPCORDtri(iptri,1)
            xi2tri=pGPCORDtri(iptri,2)
            pNNtri(iptri,1) = xi1tri
            pNNtri(iptri,2) = xi2tri
            pNNtri(iptri,3) = 1.0d0-xi1tri-xi2tri
        end do
        AB_edge = (/ (coords(1,1)-coords(1,2)) , (coords(2,1)-coords(2,2)), (coords(3,1)-coords(3,2)) /)	 
        BC_edge = (/ (coords(1,2)-coords(1,3)) , (coords(2,2)-coords(2,3)), (coords(3,2)-coords(3,3)) /)
        CD_edge = (/ (coords(1,3)-coords(1,4)) , (coords(2,3)-coords(2,4)), (coords(3,3)-coords(3,4)) /)
        AC_edge = (/ (coords(1,1)-coords(1,3)) , (coords(2,1)-coords(2,3)), (coords(3,1)-coords(3,3)) /)
        BD_edge = (/ (coords(1,2)-coords(1,4)) , (coords(2,2)-coords(2,4)), (coords(3,2)-coords(3,4)) /)
    
        N1 = cross(AB_edge,BC_edge) ! Gives normal of the face defined by nodes (/1(A),2(B),3(C))
        N1 = N1/norm(N1)
        N2 = cross(BD_edge,CD_edge) ! Gives normal of the face defined by nodes (/2(B),3(C),4(D))
        N2 = N2/norm(N2)
        N3 = cross(AC_edge,CD_edge) ! Gives normal of the face defined by nodes (/1(A),3(C),4(D))
        N3 = N3/norm(N3)
        N4 = cross(AB_edge,BD_edge) ! Gives normal of the face defined by nodes (/1(A),2(B),4(D))
        N4 = N4/norm(N4)

        if (ABS(N1(3)).eq.1) then
            NODES = (/1, 2, 3/)
            DetjTri  =(coords(1,1)-coords(1,3))*(coords(2,2)-coords(2,3))& 
                -(coords(2,1) -coords(2,3))*(coords(1,2)-coords(1,3))
            if (DetjTri<0) then
                NODES = (/ 3,2,1/)
            end if 
    !					write(*,*) "element with N1 normal along x-direction: ", jElem
    !					write(*,*) "N1 normal : ", N1
    !					write(*,*) "Jacobian : ", DetjTri
        elseif (ABS(N2(3)).eq.1) then
            NODES = (/2, 3, 4/)
            DetjTri  =(coords(1,2)-coords(1,4))*(coords(2,3)-coords(2,4))& 
                -(coords(2,2) -coords(2,4))*(coords(1,3)-coords(1,4))
            if (DetjTri<0) then
                NODES = (/ 4,3,2/)
            end if 
    !					write(*,*) "element with N2 normal along x-direction: ", jElem
    !					write(*,*) "N2 normal : ", N2
    !					write(*,*) "Jacobian : ", DetjTri
    !					write(*,*) "NODES : ", NODES
        elseif (ABS(N3(3)).eq.1) then
            NODES = (/1, 3, 4/)
            DetjTri  =(coords(1,1)-coords(1,4))*(coords(2,3)-coords(2,4))& 
                -(coords(2,1) -coords(2,4))*(coords(1,3)-coords(1,4))
            if (DetjTri<0) then
                NODES = (/ 4,3,1/)
            end if 
    !					write(*,*) "element with N3 normal along x-direction: ", jElem
    !					write(*,*) "N3 normal : ", N3
    !					write(*,*) "Jacobian : ", DetjTri
    !					write(*,*) "NODES : ", NODES
        elseif (ABS(N4(3)).eq.1) then
            NODES = (/1, 2, 4/)
            DetjTri  =(coords(1,1)-coords(1,4))*(coords(2,2)-coords(2,4))& 
                -(coords(2,1) -coords(2,4))*(coords(1,2)-coords(1,4))
            if (DetjTri<0) then
                NODES = (/ 4,2,1/)
            end if 
    !					write(*,*) "element with N4 normal along x-direction: ", jElem
    !					write(*,*) "N4 normal : ", N4
    !					write(*,*) "Jacobian : ", DetjTri
        end if
! -------------------------------------------- Application of flux vector -----------------------------------------------------
!					write(*,*) "RHS", RHS
        if (ANY(jElem.eq.FrontEle)) then
!            DO iptri=1,iGPtri
!                DO ni=1,size(NODES)
!                    nj = NODES(ni)
!                    RHS(nj) = RHS(nj) - pWTtri(iptri)*ABS(DetjTri)*kEpZero*kEpR*pNNtri(iptri,ni)*u(nj)*pa
!                    DO nk=1,size(NODES)
!                        nl = NODES(nk)
!                        AMATRX(nj,nl) = AMATRX(nj,nl) + pWTtri(iptri)*ABS(DetjTri)*kEpZero*kEpR*pNNtri(iptri,ni)*pNNtri(iptri,nk)*pa
!                    END DO
                    
!                END DO
!            END DO 		! ------------------------ iptri-loop ------------------------
        else if (ANY(jElem.eq.BackEle)) then
!            DO iptri=1,iGPtri
!                DO ni=1,size(NODES)
!                    nj = NODES(ni)
!!                    pRobin = 0.0d0
!                    RHS(nj) = RHS(nj) + pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*pRobin
!!                    DO nk=1,size(NODES)
!!                        nl = NODES(nk)
!!                        AMATRX(nj,nl) = AMATRX(nj,nl) + pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*pNNtri(iptri,nk)*pRobin
!!                    END DO
                    
!                END DO
!            END DO 		! ------------------------ iptri-loop ------------------------
        end if
    end if 		! ------------------------ Boundary jElem-loop ------------------------
!		if (jElem==56618) then	
!			write(*,*) "AMATRX", AMATRX
!		end if 
    if(jelem==26516) then
        if (kinc==1) then			
        !			write(*,*) "AMATRX", AMATRX
        !			write(*,*) "RHS:	",RHS(1:NDOFEL)
        !			write(*,*) "RHS_test:	",RHS_test(1:NDOFEL)
        !			write(*,*) "RHS_num:	",RHS_num(1:NDOFEL)
        !			write(*,*) "RHS_diff:	",RHS_num(1:NDOFEL)-RHS_test(1:NDOFEL)
        !			write(*,*)
        end if
    end if
    ! ============================================
    ! element outgoing note
    !write(*,*)" "
    !write(*,*)"UEL out - ELEMNO",jELEM
    !write(*,*)" "
    ! ============================================
!	if (JELEM == 1) then
!		write(*,*)" OUTGOING"
!           	write(*,*)"element no: ",JELEM
!            	write(*,*)" "
!	    	write(*,*)"NODAL TEMPS ", U
!	    	write(*,*)"PREDEF ", predef(2,2,:)
!	    	write(*,*)"RHS ", RHS
!	    	write(*,*)"AMATRX ", AMATRX
!	    end if

return
end subroutine UEL
!-------------------------------------------------------------------------------
!===============================================================================
