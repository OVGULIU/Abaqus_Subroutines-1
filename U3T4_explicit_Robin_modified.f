!DEC$ FREEFORM
!===============================================================================
! DEVELOPER: EDGAR HUSSER
! YEAR: 2017
!===============================================================================

! include of ABAQUS subroutines
!INCLUDE 'VEXTERNALDB.f'
!INCLUDE 'UVARM.f'
!INCLUDE 'SDVINI.f'

module functions

    implicit none
    !---------------------------------------------------------------------------          
    public :: sgn,norm,trace,det,inv,dot,ddot,dya,cross,matvec,matmat
    
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
    
    !============================================================================
    ! scalar/dot product of two vectors of the same dimension
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
   !============================================================================

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


module constitutive_relations

    implicit none
    !---------------------------------------------------------------------------          
    public :: STRESSES_CONCEN_COUPLED
    
    contains

    !===========================================================================
    ! linear elasticity (Hooke's law)
    subroutine STRESSES_CONCEN_COUPLED(S,Ee,ElecDisp,pQf,pID,pGM,pLAM,pEPSILONZERO,pEPSILONR,pEMCoup, pZ, cSat)
       
        use functions, only : trace,dya,dot
        
        implicit none
        !---------------------------------------------------------------------------  
        ! declaration
        
        double precision, dimension(:,:), intent(inout) :: S
        double precision, dimension(:,:), intent(in) :: Ee,pID
		double precision, dimension(:), intent(in) :: ElecDisp
        double precision, intent(in) :: pGM,pLAM,pEPSILONZERO,pEPSILONR,pEMCoup, pZ, cSat, pQf
        !---------------------------------------------------------------------------

        !S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID - (pEMCoup/pZ)*pQf*pID
		!write(*,*)"S",S
		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID + (1.d0/(pEPSILONZERO*pEPSILONR))*(dya(ElecDisp,ElecDisp) - 0.5d0*(dot(ElecDisp,ElecDisp))*pID)- (pEMCoup/pZ)*pQf*pID
!		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID
!		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID + (1.d0/(pEPSILONZERO*pEPSILONR))*(dya(ElecDisp,ElecDisp) - 0.5d0*(dot(ElecDisp,ElecDisp))*pID)
                !---------------------------------------------------------------------------       
            
            end subroutine STRESSES_CONCEN_COUPLED
    !===========================================================================

end module constitutive_relations
!===============================================================================
! USER SUBROUTINE - UPDATE EXTERNAL DATABASE
	SUBROUTINE VEXTERNALDB(lOp, i_Array, niArray, r_Array, nrArray)

	!include 'vaba_param.inc'
! ------Contents of i_Array------
	integer, parameter :: i_int_nTotalNodes	= 1
	integer, parameter :: i_int_nTotalElements = 2
	integer, parameter :: i_int_kStep = 3
	integer, parameter :: i_int_kInc = 4
	integer, parameter :: i_int_iStatus = 5
	integer, parameter :: i_int_lWriteRestart = 6 

! ------Possible values for lOp argument------
	integer, parameter :: j_int_StartAnalysis = 0
	integer, parameter :: j_int_StartStep = 1
	integer, parameter :: j_int_SetupIncrement = 2
	integer, parameter :: j_int_StartIncrement = 3
	integer, parameter :: j_int_EndIncrement = 4
	integer, parameter :: j_int_EndStep =5
	integer, parameter :: j_int_EndAnalysis = 6 

! ------Possible values i_Array(i_int_iStatus)------
	integer, parameter :: j_int_Continue = 0
	integer, parameter :: j_int_TerminateStep = 1
	integer, parameter :: j_int_TerminateAnalysis = 2

! ------Contents of r_Array------
	integer, parameter :: i_flt_TotalTime = 1
	integer, parameter :: i_flt_StepTime = 2
	integer, parameter :: i_flt_dtime = 3 
!
	integer, intent(in ) :: lOp,nrArray,niArray
	integer, intent(in ), dimension(niArray) :: i_Array
	double precision, intent(in ), dimension(nrArray) :: r_Array

	integer :: kstep,kInc
    logical :: I_EXIST
    character*256 :: JOBNAME
	character*256 :: filename
	kstep = i_Array(i_int_kStep)
	kInc = i_Array(i_int_kInc)
    
! ------ START OF THE ANALYSIS ------
	if (lOp .eq. j_int_StartAnalysis) then
    call VGETJOBNAME(JOBNAME,LENJOBNAME)
    filename = '/home/cerecam/Desktop/LimitReached'// trim(JOBNAME) // '.inp'
        INQUIRE(FILE=filename,EXIST=I_EXIST)
        if (I_EXIST) then
            open(unit=107, file=filename)            
            close(UNIT=107,STATUS='DELETE')
            write(*,*) " !!! LimitReached"// trim(JOBNAME) // ".inp Deleted !!"
        end if
! ------ Start of the step ------
	else if (lOp .eq. j_int_StartStep) then
! ------ Setup the increment ------
	else if (lOp .eq. j_int_SetupIncrement) then
! ------ Start of increment ------
	else if (lOp .eq. j_int_StartIncrement) then
! ------ End of increment ------
	else if (lOp .eq. j_int_EndIncrement) then
! ------ End of the step ------
	else if (lOp .eq. j_int_EndStep) then
	else if (lOp .eq. j_int_EndAnalysis) then
	end if

	return
	end subroutine VEXTERNALDB



!===============================================================================
!-------------------------------------------------------------------------------
! USER SUBROUTINE - USER ELEMENT

! abaqus interface VUEL
SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars, &
                energy, &
                nnode,ndofel,props,nprops,jprops,njprops, &
                coords,ncrd,u,du,v,a, &
                jtype,jElem, &
                time,period,dtimeCur,dtimePrev,kstep,kinc, &
                lflags, &
                dMassScaleFactor, &
                predef,npredef, &
                jdltyp, adlmag)
                
    use functions
    use constitutive_relations
    
    implicit none ! use explicit=both and output_precision=full
    !include 'vaba_param.inc'
    ! if implicit none is not used and 'vaba_param.inc' is included - implicit declaration of variables:
    ! all variables for letters a-h, o-z is real and i-n are integers (i,j,k,l,m,n).
    !-------------------------------------------------------------------------------   
    ! declaration
    
    ! parameters - numbers
    !parameter ( zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, four=4.d0, six=6.d0)
    double precision, parameter :: zero=0.d0
    double precision, parameter :: half=0.5d0
    double precision, parameter :: one=1.d0
    double precision, parameter :: two=2.d0
    double precision, parameter :: three=3.0d0
    double precision, parameter :: four=4.d0
    double precision, parameter :: six=6.d0
    double precision, parameter :: factorStable=0.99d0
    double precision, parameter :: pi=3.1415926535897932
    
    ! parameters - problem specification
    !parameter ( iGP=1, iCORD=3, iNODE=4 )
    integer, parameter :: iGP=1		!Number of Gauss Points
    integer, parameter :: iCORD=3	!Degrees of freedom (mechanical)
    integer, parameter :: ICORDTOTAL=4  !Degrees of freedom (total per node (3 mech; 1 temp))
    integer, parameter :: iNODE=4	!Number of nodes per element
    integer, parameter :: iGPtri=3          !Number of Gauss Points for triangle surface elements
    
    ! Abaqus variables
    ! ================  
    
    ! predefined parameters - operational code keys
    integer, parameter :: jMassCalc = 1
    integer, parameter :: jIntForceAndDtStable = 2
    integer, parameter :: jExternForce = 3
    
    ! predefined parameters - flag indices
    integer, parameter :: iProcedure = 1
    integer, parameter :: iNlgeom = 2
    integer, parameter :: iOpCode = 3
    integer, parameter :: nFlags = 3
    
    !  predefined parameters - procedure flags
    integer, parameter :: jDynExplicit = 17
    
    ! predefined parameters - energy array indices
    integer, parameter :: iElPd = 1
    integer, parameter :: iElCd = 2
    integer, parameter :: iElIe = 3
    integer, parameter :: iElTs = 4
    integer, parameter :: iElDd = 5
    integer, parameter :: iElBv = 6
    integer, parameter :: iElDe = 7
    integer, parameter :: iElHe = 8
    integer, parameter :: iElKe = 9
    integer, parameter :: iElTh = 10
    integer, parameter :: iElDmd = 11
    integer, parameter :: iElDc = 12
    integer, parameter :: nElEnergy = 12

    ! predefined parameters - time indices
    integer, parameter :: iStepTime = 1
    integer, parameter :: iTotalTime = 2
    integer, parameter :: nTime = 2
    
    ! predefined parameters - predefined variables indices
    integer, parameter :: iPredValueNew = 1
    integer, parameter :: iPredValueOld = 2
    integer, parameter :: nPred = 2

    ! variables passed in for information
    integer, intent(in ) :: nblock                              ! number of user elements to be processed in this call to VUEL
    double precision, intent(in ) :: dtimeCur                   ! current time increment
    double precision, intent(in ) :: dtimePrev                  ! previous time increment
    double precision, intent(in ) :: period                     ! time period of the current step
    integer, intent(in ) :: ndofel                              ! number of dofs in the element
    integer, intent(in ) :: nsvars                    	        ! user-defined number of solution-dependent variables            
    integer, intent(in ) :: nprops                              ! user-defined number of real property values
    integer, intent(in ) :: njprops                             ! user-defined number of integer property values 
    integer, intent(in ) :: ncrd                                ! max. number of user-defined coordinates at any nodal point
    integer, intent(in ) :: nnode                               ! user-defined number of nodes on the element
    integer, intent(in ) :: jtype                               ! integer defining the element type
    integer, intent(in ) :: kstep 	                            ! current step number
    integer, intent(in ) :: kinc 	                            ! current increment number
    integer, intent(in ) :: npredef                             ! number of predefined field variables (incl. temperature)
    integer, intent(in ) :: jdltyp                              ! specifies the load type
    
    double precision, intent(in ), dimension(nprops) :: props   ! real property values
    integer, intent(in ), dimension(njprops) :: jprops          ! integer property values
    double precision, dimension(nblock,nnode,ncrd) :: coords    ! block of original coordinates of nodes (reference configuration)
    integer, intent(in ), dimension(nblock) :: jElem 	        ! block of element numbers
    double precision, intent(in ), dimension(nblock) :: adlmag 	! block of total load magnitute of type jdltyp
    double precision, intent(in ), dimension(nblock,nnode,npredef,nPred) :: predef 	! block of values of predefined field variables (only uncoupled)
    integer, intent(in ), dimension(nFlags) :: lflags 	        ! array containing the flags of current solution procedure
    double precision, intent(in ), dimension(nblock) :: dMassScaleFactor            ! block containing mass scale factors (for each element)
    double precision, intent(in ), dimension(nTime) :: time     ! current value of step time
    double precision, intent(in ), dimension(nblock,ndofel) :: u    ! block of nodal solution variables (displacement, temperature, etc.)
    double precision, intent(in ), dimension(nblock,ndofel) :: du   ! block of incremental values (displacement, temperature, etc.)
    double precision, intent(in ), dimension(nblock,ndofel) :: v    ! block of time rates of nodal solution variables
    double precision, intent(in ), dimension(nblock,ndofel) :: a    ! block of acclerations of nodal solution variables
    
    ! variables to be defined (block-wise)
    double precision, intent(inout), dimension(nblock,ndofel) :: rhs 	        ! element contribution to right-hand-side vector
    double precision, intent(inout), dimension(nblock,ndofel,ndofel) :: amass  ! element contribution to mass matrix (symmetric)
    double precision, intent(inout), dimension(nblock) :: dtimeStable  	    ! scalar value, upper limit of time increment (Courant condition)
    double precision, intent(inout), dimension(nblock,nsvars) :: svars  	    ! solution-depending state variables
    double precision, intent(inout), dimension(nblock,nElEnergy) :: energy		! element energy quantities 
    
    ! user variables
    ! ==================
    
    double precision :: pGPCORD(iGP,iCORD)  ! integration point coordinates
    double precision :: pQUAD               ! factor for quadrature rule
    double precision :: pWT(iGP)            ! integration point weights
    double precision :: pNN(iGP,iNODE)      ! shape function matrix (interpolation matrix)
    double precision :: pID(iCORD,iCORD)    ! identity tensor

    ! user variables- For triangular elements
    ! ==================
    
    double precision :: pGPCORDtri(iGPtri,iCORD-1)  ! integration point coordinates
    double precision :: pQUADtri               ! factor for quadrature rule
    double precision :: pWTtri(iGPtri)            ! integration point weights
    double precision :: pNNtri(iGPtri,iNODE-1)      ! shape function matrix (interpolation matrix)
    
    double precision :: amass_row_sum       ! sum of array row
    double precision :: cd,pd,pd_min
    double precision :: Xp(iCORD),Xa(iCORD),Xb(iCORD),Xc(iCORD)
    double precision :: Xba(iCORD),Xca(iCORD),Xpa(iCORD)
    double precision :: pNb(iCORD),pN(iCORD)

    ! include others
    !INCLUDE 'COMMON.f'

    !===============================================================================

    double precision :: xi1,xi2,xi3         ! natural coordinates
    double precision :: xi1tri,xi2tri,xi3tri         ! natural coordinates
    double precision :: X1(iNODE)           ! physical coordinates
    double precision :: X2(iNODE)           ! physical coordinates
    double precision :: X3(iNODE)           ! physical coordinates
    double precision :: dNdXi1(iNODE)       ! shape function derivatives
    double precision :: dNdXi2(iNODE)       ! shape function derivatives
    double precision :: dNdXi3(iNODE)       ! shape function derivatives
    double precision :: dNdX1(iGP,iNODE)    ! shape function derivatives
    double precision :: dNdX2(iGP,iNODE)    ! shape function derivatives
    double precision :: dNdX3(iGP,iNODE)    ! shape function derivatives
    double precision :: dX1dxi1             ! derivatives
    double precision :: dX1dxi2             ! derivatives
    double precision :: dX1dxi3             ! derivatives
    double precision :: dX2dxi1             ! derivatives
    double precision :: dX2dxi2             ! derivatives
    double precision :: dX2dxi3             ! derivatives
    double precision :: dX3dxi1             ! derivatives
    double precision :: dX3dxi2             ! derivatives
    double precision :: dX3dxi3             ! derivatives
    double precision :: JJ(iCORD,iCORD)     ! Jacobi matrix
    double precision :: detJ(iGP)           ! Jacobi-determinant (reference)
    double precision :: DetjTri           ! Jacobi-determinant triangular faces (reference)
    double precision :: Uarray(iCORD)       ! nodal displacement in each dimension
    double precision :: CoNODE(iNODE)       ! nodal concentration 1.4498201477996281E-007s
    double precision :: Co           ! Concentration inetrpolated in element
    double precision :: pCo          ! Concentration mobile interpolated in element
    double precision :: pQf       ! Free charge density
    double precision :: phi          ! nodal Electrical potentials
    double precision :: gCo(iCORD)          ! Grad of concentrations
    double precision :: Ux(iNODE)           ! nodal displacement in x
    double precision :: Uy(iNODE)           ! nodal displacement in y
    double precision :: Uz(iNODE)           ! nodal displacement in z
    double precision :: H(iCORD,iCORD)      ! displacement gradient
    double precision :: F(iCORD,iCORD)      ! deformation gradient
    double precision :: Ee(iCORD,iCORD)     ! small /Gree-Lagrange strain tensor
    double precision :: SIG(iCORD,iCORD)    ! Cauchy stress tensor
    double precision :: S(iCORD,iCORD)      ! 2. PK
    double precision :: P(iCORD,iCORD)      ! 1. PK
    double precision :: pELECFIELD(iCORD)    ! Electric Field vector from grad(phi)
    double precision :: ELECDISP(iCORD)    ! Electric Displacement vector calculated from elec field vector
    double precision :: J(iCORD)    ! flux vector for diffusion problem
    !--------------------------Constants used in Constitutive models-------------------------------------

    double precision :: pEM		! Young's Modulus
    double precision :: pNU		! Poisson's ratio
    double precision :: pRHO            ! density
    double precision :: pLAM		! Lambda (Lame parameter 1)
    double precision :: pGM		! Mu (shear modulus, Lame parameter 2)
    double precision :: pEPSILONR	! Vacuum permittivitty
    double precision :: pEPSILONZERO	! Relative Permittivity
    double precision :: pRTHETA		! Gas constant times temp 
    double precision :: pF		! Faraday's constant 
    double precision :: pZ		! Valence mobile species 
    double precision :: pEMCoup		! Electromechanical coupling parameter
    double precision :: pDif		! Diffusivity
    double precision :: pDifnew
    double precision :: Csat		! Saturated concentration
    double precision :: pCent(iCORD)	! Centroid
    double precision :: pVolume		! Volume of element
    double precision :: pSED		! SED
    double precision :: pa1		! stability term RHS
    double precision :: pRobin		! stability term AMASS
    double precision :: Pe		! Peclet number
    double precision :: Courant		! Courant number
    double precision :: Sigma_K     ! Stability parameter [ref. Tu et al, Compt Physics Comm]
    double precision :: pA_VECTOR(iCORD)	! Stabilization vector [ref. Tu et al, Compt Physics Comm]
    double precision :: V_DASH		! Additive diffusion
    double precision :: elesize_sum	! Element size intermediate term
    double precision :: Elesize		! Element size in veocity field direction (Radius of equivalent volume circle)
    double precision :: pCoTotal	! integration within domain
    double precision :: pStability(iCORD)	! stability        

!--------------------------Constants used in modified PNP model-------------------------------------

    double precision :: pV_i		! Volume concentration of any mole of a species
    double precision :: pV_tot      ! Maximum volume concentration allowed
    double precision :: pDensVolFrac! Total density volume fraction of mobile species
    double precision :: pGradRho	! Grad of total volume
    double precision :: pNa		    ! Avogadro's constant
    double precision :: pPi		    ! Constant of pi
    double precision :: pRi		    ! Radius of concentration molecules
    double precision :: pImmobileConc		    ! Radius of concentration molecules
                
    ! integer
    integer :: kblock,ip,nn,ni,nj,i


    integer :: CordRange(iCORD) = (/(i, i=1,iCORD, 1)/)
    integer :: NODERange(iNODE) = (/(i, i=1,iNODE, 1)/)
    
    !--------------------------Constants used in Robin BC implementation -------------------------------------
    !double precision :: vMeq
    integer :: dofni(iCORD)        ! current dof
    integer :: dofniT
    integer :: dofnj(iCORD)        ! current dof
    integer :: dofnjT
        
    double precision :: AB_edge(iCORD)           ! Edge between nodes 1 and 2
    double precision :: BC_edge(iCORD)           ! Edge between nodes 2 and 3
    double precision :: CD_edge(iCORD)           ! Edge between nodes 3 and 4
    double precision :: AC_edge(iCORD)           ! Edge between nodes 1 and 3
    double precision :: BD_edge(iCORD)           ! Edge between nodes 2 and 4

    double precision :: N1(iCORD)           ! Normal described by face given by vertex: 1,2,3
    double precision :: N2(iCORD)           ! Normal described by face given by vertex: 2,3,4
    double precision :: N3(iCORD)           ! Normal described by face given by vertex: 1,3,4
    double precision :: N4(iCORD)           ! Normal described by face given by vertex: 1,2,4
    double precision :: N_all(iCORD)           ! Normals described by different faces

    double precision :: NODES(iCORD)           ! physical coordinates
    ! integer
    integer :: iptri,Filesize
    ! Allocatable arrays
!    integer, DIMENSION(:), ALLOCATABLE :: FrontEle 
!    integer, DIMENSION(:), ALLOCATABLE :: BackEle 
    integer, DIMENSION(:), ALLOCATABLE :: Z0Ele 
    integer, DIMENSION(:), ALLOCATABLE :: Z1Ele
    
    logical :: I_EXIST        
    character*256 :: JOBNAME
    character*256 :: filename
    integer :: LENJOBNAME
                
    if (jtype.eq.34) then 
    
        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetFront.inp',status='old')!
        READ(107,*) Filesize
        Allocate ( Z1Ele(Filesize) )
        close(107)

        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetFront.csv',status='old')!
        READ(107,*) Z1Ele
        close(107)

        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetBack.inp',status='old')!
        READ(107,*) Filesize
        Allocate ( Z0Ele(Filesize) )
        close(107)

        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetBack.csv',status='old')!
        READ(107,*) Z0Ele
        close(107)
        
!        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/ShortExperimental/NumberZ0_ELEMENTS.inp',status='old')!
!        READ(107,*) Filesize
!        Allocate ( Z0Ele(Filesize) )
!        close(107)

!        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/ShortExperimental/Z0_ELEMENTS.csv',status='old')!
!        READ(107,*) Z0Ele
!        close(107)
        
!        Filesize=655
!        Allocate ( Z1Ele(Filesize) )

!        open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short1/Z1_ELEMENTS.csv',status='old')!
!        READ(107,*) Z1Ele
!        close(107)
    
        pEM  = props(1)
        pNU  = props(2)
        pRHO = props(3)
        pEPSILONR  = props(4)
        pEPSILONZERO  = props(5)
        pRTHETA = props(6)
		pF = props(7)
		pZ = props(8)
		pEMCoup = props(9)
		pDif = props(10)
		cSat = props(11)
		pa1 = props(12)
		pRobin = props(13)
        
        pGM  = half*pEM/(one+pNU)
        pLAM = (pEM*pNU)/((one+pNU)*(one-two*pNU))

!--------------------------Parameters used in modified PNP model-------------------------------------
        pNa = 602.2140857
        pPi = 3.14159265358979311
        pRi = 0.502
        !pImmobileConc = 1.8e-3
        !pV_tot = 4.0d0/3.0d0*pPi*(pRi**3)*pNa*(pImmobileConc)
        pV_tot = 4.0
        pImmobileConc = pV_tot / (4.0d0/3.0d0*pPi*(pRi**3)*pNa) 
           
        ! integration point coordinates and weights --------------------------------
        if (iGP==1) then ! HUGHES - The Finite Element Method 1987 (p. 174)

            pGPCORD(1,:) = (/ one/four, one/four, one/four /)
                
            pWT(1) = one
                
            pQUAD = one/six
        end if

        !---------------------------------------------------------------------------
            
        ! shape function for each ip -----------------------------------------------
        do ip=1,iGP

            xi1=pGPCORD(ip,1)
            xi2=pGPCORD(ip,2)
            xi3=pGPCORD(ip,3)

            if (iNODE==4) then ! cf. WRIGGERS - Nonlinear Finite Elemente Methods 2008 (p. 120)

                pNN(ip,1) = one-xi1-xi2-xi3
                pNN(ip,2) = xi1
                pNN(ip,3) = xi2
                pNN(ip,4) = xi3

            else
                stop "Error in computation of shape functions. The number of nodes does not conform with the element type (4 node tetrahedral element)."
            end if

        end do
        !---------------------------------------------------------------------------
            
        ! identity matrix ----------------------------------------------------------
        pID=zero
        forall(i=1:iCORD) pID(i,i)=one
        !---------------------------------------------------------------------------
        
        if (lflags(iOpCode).eq.jMassCalc) then ! compute mass matrix
            
            ! loop over element block
            do kblock=1,nblock ! ---------------------------------------------------
            
                ! loop over all integration points (computation of FE variables)
                do ip=1,iGP ! ------------------------------------------------------
        
                    ! get solution-dependent state variables (history variables)
        
                    ! natural coordinates of current ip
                    xi1 = pGPCORD(ip,1)
                    xi2 = pGPCORD(ip,2)
                    xi3 = pGPCORD(ip,3)
        
                    ! coordinate vectors of current element
                    X1 = coords(kblock,:,1)
                    X2 = coords(kblock,:,2)
                    x3 = coords(kblock,:,3)
    
                    ! --------------------------------------------------------------
                    if (iNODE==4) then
    
                        ! derivatives of shape functions with respect to natural coordinates                
                        dNdXi1(1) = -one
                        dNdXi2(1) = -one
                        dNdXi3(1) = -one
        
                        dNdXi1(2) =  one
                        dNdXi2(2) =  zero
                        dNdXi3(2) =  zero
        
                        dNdXi1(3) =  zero
                        dNdXi2(3) =  one
                        dNdXi3(3) =  zero
        
                        dNdXi1(4) =  zero
                        dNdXi2(4) =  zero
                        dNdXi3(4) =  one
                     
                    else 
                        stop "Error in computation of shape function derivatives. The number of nodes does not conform with the element type (4 node tetrahedral element)."     
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
                    do nn=1,iNODE
                        dNdX1(ip,nn) = one/detJ(ip)*( (dX2dxi2*dX3dxi3-dX3dxi2*dX2dxi3)*dNdXi1(nn) &
                                                    + (dX3dxi1*dX2dxi3-dX2dxi1*dX3dxi3)*dNdXi2(nn) &
                                                    + (dX2dxi1*dX3dxi2-dX3dxi1*dX2dxi2)*dNdXi3(nn) )
                        dNdX2(ip,nn) = one/detJ(ip)*( (dX3dxi2*dX1dxi3-dX1dxi2*dX3dxi3)*dNdXi1(nn) &
                                                    + (dX1dxi1*dX3dxi3-dX3dxi1*dX1dxi3)*dNdXi2(nn) &
                                                    + (dX3dxi1*dX1dxi2-dX1dxi1*dX3dxi2)*dNdXi3(nn) )
                        dNdX3(ip,nn) = one/detJ(ip)*( (dX1dxi2*dX2dxi3-dX2dxi2*dX1dxi3)*dNdXi1(nn) &
                                                    + (dX2dxi1*dX1dxi3-dX1dxi1*dX2dxi3)*dNdXi2(nn) &
                                                    + (dX1dxi1*dX2dxi2-dX2dxi1*dX1dxi2)*dNdXi3(nn) )
                    end do
    
                end do ! -----------------------------------------------------------
                
                !amass(kblock,1:ndofel,1:ndofel)=zero
                ! loop over all integration points (computation of mass matrix)
                do ip=1,iGP ! ------------------------------------------------------

                    ! summation over node_i
                    do ni=1,iNODE !-----------------------------loop-i--------------
        		
                        ! current node dof
                        dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
						dofniT = iCORDTOTAL*ni
                        ! summation over node_i
                        do nj=1,iNODE !-------------------------loop-j--------------
            
                            ! current node dof
                            dofnj(1:iCORD) = 1+iCORDTOTAL*((nj*1)-1)+(CordRange-1)
                            dofnjT = iCORDTOTAL*nj
    
                            ! regular mass matrix
                            amass(kblock,dofni,dofnj) = amass(kblock,dofni,dofnj) &
                                + pQUAD*pWT(ip)*detJ(ip)*pRHO*matmat(pNN(ip,ni)*Pid,pNN(ip,nj)*Pid)
                            ! Capacitence matrix calculation
                            amass(kblock,dofniT,dofnjT) = amass(kblock,dofniT,dofnjT) &
                                + pQUAD*pWT(ip)*detJ(ip)*pNN(ip,ni)*pNN(ip,nj)
                        end do !--------------------------end-loop-j----------------
    
                    end do !------------------------------end-loop-i----------------
    
                end do ! ----------------------------end-loop-ip--------------------
                
                ! mass lumbing
                do i=1,ndofel
                    amass_row_sum = sum(amass(kblock,i,:))
                    amass(kblock,i,:) = zero
                    amass(kblock,i,i) = amass_row_sum
                end do
               

            end do !----------------------------------------------------------------
                
        elseif (lflags(iOpCode).eq.jIntForceAndDtStable) then !compute internal force + stable time increment

            ! loop over element block
            do kblock=1,nblock ! ---------------------------------------------------

    
                ! loop over all integration points (computation of FE variables)
                do ip=1,iGP ! ------------------------------------------------------
        
                    ! get solution-dependent state variables (history variables)
        
                    ! natural coordinates of current ip
                    xi1 = pGPCORD(ip,1)
                    xi2 = pGPCORD(ip,2)
                    xi3 = pGPCORD(ip,3)
        
                    ! coordinate vectors of current element
                    X1 = coords(kblock,:,1)
                    X2 = coords(kblock,:,2)
                    x3 = coords(kblock,:,3)
    
                    ! --------------------------------------------------------------
                    if (iNODE==4) then
    
                        ! derivatives of shape functions with respect to natural coordinates                
                        dNdXi1(1) = -one
                        dNdXi2(1) = -one
                        dNdXi3(1) = -one
        
                        dNdXi1(2) =  one
                        dNdXi2(2) =  zero
                        dNdXi3(2) =  zero
        
                        dNdXi1(3) =  zero
                        dNdXi2(3) =  one
                        dNdXi3(3) =  zero
        
                        dNdXi1(4) =  zero
                        dNdXi2(4) =  zero
                        dNdXi3(4) =  one
                     
                    else 
                        stop "Error in computation of shape function derivatives. The number of nodes does not conform with the element type (4 node tetrahedral element)."     
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
                    do nn=1,iNODE
                        dNdX1(ip,nn) = one/detJ(ip)*( (dX2dxi2*dX3dxi3-dX3dxi2*dX2dxi3)*dNdXi1(nn) &
                                                    + (dX3dxi1*dX2dxi3-dX2dxi1*dX3dxi3)*dNdXi2(nn) &
                                                    + (dX2dxi1*dX3dxi2-dX3dxi1*dX2dxi2)*dNdXi3(nn) )
                        dNdX2(ip,nn) = one/detJ(ip)*( (dX3dxi2*dX1dxi3-dX1dxi2*dX3dxi3)*dNdXi1(nn) &
                                                    + (dX1dxi1*dX3dxi3-dX3dxi1*dX1dxi3)*dNdXi2(nn) &
                                                    + (dX3dxi1*dX1dxi2-dX1dxi1*dX3dxi2)*dNdXi3(nn) )
                        dNdX3(ip,nn) = one/detJ(ip)*( (dX1dxi2*dX2dxi3-dX2dxi2*dX1dxi3)*dNdXi1(nn) &
                                                    + (dX2dxi1*dX1dxi3-dX1dxi1*dX2dxi3)*dNdXi2(nn) &
                                                    + (dX1dxi1*dX2dxi2-dX2dxi1*dX1dxi2)*dNdXi3(nn) )
                    end do
    
                end do ! -----------------------------------------------------------       
    
                rhs(kblock,1:ndofel)=zero
                energy(kblock,iElIe)=zero
                ! loop over all integration points (computation of residuum)
                do ip=1,iGP ! ------------------------------------------------------
    
                    ! displacement gradient
                    Ux = u(kblock,1:iCORD*iNODE:iCORD)
                    Uy = u(kblock,2:iCORD*iNODE:iCORD)
                    Uz = u(kblock,3:iCORD*iNODE:iCORD)
                    H = zero
                    gCo = zero
                    pELECFIELD = zero
                    do ni=1,iNODE
                        dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
                        Uarray(1:iCORD)=u(kblock, dofni)
                        H = H + dya(Uarray,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
                        CoNODE(ni) = u(kblock, (iCORDTOTAL*ni))
                        gCo = gCo + ( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )*CoNODE(ni)
                        !Electric Field
                        pELECFIELD = pELECFIELD - ((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))*predef(kblock,ni,2,1)
                    end do
                    pCo = dot(pNN(ip,:),CoNODE)
                    if (pCo.lt.0.0d0) then
                        pCo = 0.0d0
                    end if
    				pV_i = 4.0d0/3.0d0*pPi*(pRi**3)*pNa*pCo
                
                    pDensVolFrac = pV_i/pV_tot
                    if (pDensVolFrac>1.0) then
!                        write(*,*) "rho exceeds to 1.0"
                    end if
                    if (pDensVolFrac<0.0) then
                        pDensVolFrac=0.0
!                        write(*,*) "rho exceeds to 0.0"
                    end if
                
                    ! small strain tensor
                    Ee = half*(transpose(H) + H)
                    
                    ElecDisp = pEPSILONZERO*pEPSILONR*pELECFIELD
					pQf = pF*((pZ*pCo)+(cSat*(1.d0)))
                    ! stress
                    call stresses_concen_coupled(S,Ee,ElecDisp,pQf,pID,pGM,pLAM,pEPSILONZERO,pEPSILONR,pEMCoup, pZ, cSat)

                    energy(kblock,iElIe)= energy(kblock,iElIe) + (detJ(ip)/6.0d0)*( ( pGM*ddot(Ee,Ee)+half*pLAM*trace(Ee)*trace(Ee)) )

                    
                    Elesize = (((3.0d0/(pi*4.0d0))*(detJ(1)/6.0d0))**(one/3.0d0))*2
                    if (ALL(pELECFIELD.eq.0.0d0)) then
                        Pe = 0.0d0
                    else
                        Pe = NORM((pZ*pF/pRTHETA*pELECFIELD*(Elesize)))/2
                    end if
                    
                    if (Pe>1.0d0) then
                        Pe = 1.0d0
                    end if
                    
                    pA_Vector = pDif*pZ*pF/pRTHETA*pELECFIELD
                    
                    if (ALL(pELECFIELD.eq.0.0d0)) then
                        sigma_k = 0.0d0
                        Courant = 0.0d0
                    else
                        sigma_k = (Elesize/(2*NORM(pA_Vector)))*Pe
                        Courant = NORM((pDif*pZ*pF/pRTHETA*pELECFIELD))*dtimeCur/Elesize
                    end if                    

    
                    ! summation over node_i

                    do ni=1,iNODE !-----------------------------loop-i--------------
                
                        ! current node dof
						call STRESSES_CONCEN_COUPLED(S,Ee,ElecDisp,pQf,pID,pGM,pLAM,pEPSILONZERO,pEPSILONR,pEMCoup, pZ, cSat)
						dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
						dofniT = iCORDTOTAL*ni
        
                        !--------------------------------------Displacement RHS--------------------------------------
						rhs(kblock,dofni) = rhs(kblock,dofni) 	+ pQUAD*pWT(ip)*detJ(ip)*(matvec(S,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))) 
!											- pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),((pEMCoup/pZ)*pQf*pID))
!                        !--------------------------------------Concentration RHS--------------------------------------
                        rhs(kblock,dofniT) = rhs(kblock,dofniT) &
                            - (pDif)*pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(-gCo)) &
                            - (pDif)*pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(((pF*pZ)/(pRTHETA)*pNN(ip,ni)*pCo*pELECFIELD))) &
                            + (pDif)*(pF*pZ)/(pRTHETA)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,(sigma_k*pA_Vector)*pa1)))
                    
!                        !-------------------------- RHS modified PNP model-------------------------------------
!                        if (pDensVolFrac>0.1) then
!                            CALL VGETJOBNAME( JOBNAME,LENJOBNAME)
!                            filename = '/home/cerecam/Desktop/LimitReached' // trim(JOBNAME) // '.inp'
!                            INQUIRE(FILE= filename ,EXIST=I_EXIST)
!                            if (.NOT. I_Exist) then                                
!                                WRITE(*,*) "Limit has been reached in ",jelem(kblock)," at concentration of ", pCo,"kinc: ", kInc
!                                open(unit=107, file=filename)
!                                WRITE(107,*) "Limit has been reached in ",jelem(kblock)," at concentration of ", pCo," kinc: ", kInc, " for job: ", trim(JOBNAME)
!                                close(107)
!                            end if
!                                    rhs(kblock,dofniT) = rhs(kblock,dofniT) &
!                            - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(-(one-pDensVolFrac)*gCo)) &
!                            - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(pDif*pNN(ip,ni)*pCo*(1/(pImmobileConc))*gCo)) &
!                            + (pF*pZ)/(pRTHETA)*(one-pDensVolFrac)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,(sigma_k*pA_Vector)*pa1)))
!                        else
!                                    rhs(kblock,dofniT) = rhs(kblock,dofniT) &
!                            - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(-(one-pDensVolFrac)*gCo)) &
!                            - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(((pF*pZ)/(pRTHETA)*pNN(ip,ni)*pCo*(one-pDensVolFrac)*pELECFIELD))) &
!                            - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(pDif*pNN(ip,ni)*pCo*(1/(pImmobileConc))*gCo)) &
!                            + (pF*pZ)/(pRTHETA)*(one-pDensVolFrac)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,(sigma_k*pA_Vector)*pa1)))
                        
!                        end if
    
                    end do !------------------------------end-loop-i----------------
    
                end do ! ----------------------------end-loop-ip--------------------

!                if (ANY(jElem(kblock).eq.FrontEle) .OR. ANY(jElem(kblock).eq.BackEle)) then
                if (ANY(jElem(kblock).eq.Z0Ele)) then            
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
                    AB_edge = (/ (coords(kblock,1,1)-coords(kblock,2,1)) , (coords(kblock,1,2)-coords(kblock,2,2)), (coords(kblock,1,3)-coords(kblock,2,3)) /)	 
                    BC_edge = (/ (coords(kblock,2,1)-coords(kblock,3,1)) , (coords(kblock,2,2)-coords(kblock,3,2)), (coords(kblock,2,3)-coords(kblock,3,3)) /)
                    CD_edge = (/ (coords(kblock,3,1)-coords(kblock,4,1)) , (coords(kblock,3,2)-coords(kblock,4,2)), (coords(kblock,3,3)-coords(kblock,4,3)) /)
                    AC_edge = (/ (coords(kblock,1,1)-coords(kblock,3,1)) , (coords(kblock,1,2)-coords(kblock,3,2)), (coords(kblock,1,3)-coords(kblock,3,3)) /)
                    BD_edge = (/ (coords(kblock,2,1)-coords(kblock,4,1)) , (coords(kblock,2,2)-coords(kblock,4,2)), (coords(kblock,2,3)-coords(kblock,4,3)) /)
        
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
                        DetjTri  =(coords(kblock,1,1)-coords(kblock,3,1))*(coords(kblock,2,2)-coords(kblock,3,2))& 
                            -(coords(kblock,1,2) -coords(kblock,3,2))*(coords(kblock,2,1)-coords(kblock,3,1))
                        if (DetjTri<0) then
                            NODES = (/ 3,2,1/)
                        end if 
    !					write(*,*) "element with N1 normal along x-direction: ", jElem(kblock)
    !					write(*,*) "N1 normal : ", N1
    !					write(*,*) "Jacobian : ", DetjTri
                    elseif (ABS(N2(3)).eq.1) then
                        NODES = (/2, 3, 4/)
                        DetjTri  =(coords(kblock,2,1)-coords(kblock,4,1))*(coords(kblock,3,2)-coords(kblock,4,2))& 
                            -(coords(kblock,2,2) -coords(kblock,4,2))*(coords(kblock,3,1)-coords(kblock,4,1))
                        if (DetjTri<0) then
                            NODES = (/ 4,3,2/)
                        end if 
    !					write(*,*) "element with N2 normal along x-direction: ", jElem(kblock)
    !					write(*,*) "N2 normal : ", N2
    !					write(*,*) "Jacobian : ", DetjTri
                    elseif (ABS(N3(3)).eq.1) then
                        NODES = (/1, 3, 4/)
                        DetjTri  =(coords(kblock,1,1)-coords(kblock,4,1))*(coords(kblock,3,2)-coords(kblock,4,2))& 
                            -(coords(kblock,1,2) -coords(kblock,4,2))*(coords(kblock,3,1)-coords(kblock,4,1))
                        if (DetjTri<0) then
                            NODES = (/ 4,3,1/)
                        end if 
    !					write(*,*) "element with N3 normal along x-direction: ", jElem(kblock)
    !					write(*,*) "N3 normal : ", N3
    !					write(*,*) "Jacobian : ", DetjTri
                    elseif (ABS(N4(3)).eq.1) then
                        NODES = (/1, 2, 4/)
                        DetjTri  =(coords(kblock,1,1)-coords(kblock,4,1))*(coords(kblock,2,2)-coords(kblock,4,2))& 
                            -(coords(kblock,1,2) -coords(kblock,4,2))*(coords(kblock,2,1)-coords(kblock,4,1))
                        if (DetjTri<0) then
                            NODES = (/ 4,2,1/)
                        end if 
    !					write(*,*) "element with N4 normal along x-direction: ", jElem(kblock)
    !					write(*,*) "N4 normal : ", N4
    !					write(*,*) "Jacobian : ", DetjTri
                    end if
    ! ------------------------------------ Application of flux vector -----------------------------------------------
                    DO iptri=1,iGPtri
                        do ni=1,size(NODES)
                            nj = NODES(ni)                            
                            dofniT = iCORDTOTAL*nj
!                            if (ANY(jElem(kblock).eq.FrontEle)) then
                            if (ANY(jElem(kblock).eq.Z0Ele)) then
!                                if (pDensVolFrac>0.1) then
!                                    RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*0.0d0
!                                else                                
                                    RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*(one-pDensVolFrac)*1.5E-4
!                                end if                           
                            elseif (ANY(jElem(kblock).eq.Z1Ele)) then                               
                                RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*pDensVolFrac*(-1.0d0)*1.5E-4
                            end if
!                            if (ANY(jElem(kblock).eq.FrontEle)) then
                                
!                                RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*(-1.0d0)*3.0E-4
!                            elseif (ANY(jElem(kblock).eq.BackEle)) then
!                                if (pDensVolFrac>0.015) then
!                                    RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*0.0d0
!                                else                                
!                                    RHS(kblock,dofniT) = RHS(kblock,dofniT) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*(one-pDensVolFrac)*3.0E-4
!                                end if
!                            end if
                                    
                        END DO ! ------------------------ ni-loop ------------------------
                    END DO ! ------------------------ iptri-loop ------------------------
                end if ! ------------------------ Boundary jElem-loop ------------------------
        
                cd = sqrt( (pEM*(one-pNU))/(pRHO*(one+pNU)*(one-two*pNU)) )
                            do i=1,iNODE
!                if (i==1) then
!                    Xp = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
!                    Xa = (/coords(kblock,2,1)+u(kblock,5),coords(kblock,2,2)+u(kblock,6),coords(kblock,2,3)+u(kblock,7)/)
!                    Xb = (/coords(kblock,3,1)+u(kblock,9),coords(kblock,3,2)+u(kblock,10),coords(kblock,3,3)+u(kblock,11)/)
!                    Xc = (/coords(kblock,4,1)+u(kblock,13),coords(kblock,4,2)+u(kblock,14),coords(kblock,4,3)+u(kblock,15)/)
!                elseif (i==2) then
!                    Xp = (/coords(kblock,2,1)+u(kblock,5),coords(kblock,2,2)+u(kblock,6),coords(kblock,2,3)+u(kblock,7)/)
!                    Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
!                    Xb = (/coords(kblock,3,1)+u(kblock,9),coords(kblock,3,2)+u(kblock,10),coords(kblock,3,3)+u(kblock,11)/)
!                    Xc = (/coords(kblock,4,1)+u(kblock,13),coords(kblock,4,2)+u(kblock,14),coords(kblock,4,3)+u(kblock,15)/)
!                elseif (i==3) then
!                    Xp = (/coords(kblock,3,1)+u(kblock,9),coords(kblock,3,2)+u(kblock,10),coords(kblock,3,3)+u(kblock,11)/)
!                    Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
!                    Xb = (/coords(kblock,2,1)+u(kblock,5),coords(kblock,2,2)+u(kblock,6),coords(kblock,2,3)+u(kblock,7)/)

!                    Xc = (/coords(kblock,4,1)+u(kblock,13),coords(kblock,4,2)+u(kblock,14),coords(kblock,4,3)+u(kblock,15)/)
!                elseif (i==4) then
!                    Xp = (/coords(kblock,4,1)+u(kblock,13),coords(kblock,4,2)+u(kblock,14),coords(kblock,4,3)+u(kblock,15)/)
!                    Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
!                    Xb = (/coords(kblock,2,1)+u(kblock,5),coords(kblock,2,2)+u(kblock,6),coords(kblock,2,3)+u(kblock,7)/)
!                    Xc = (/coords(kblock,3,1)+u(kblock,9),coords(kblock,3,2)+u(kblock,10),coords(kblock,3,3)+u(kblock,11)/)
!                end if	   
                if (i==1) then
                    Xp = (/coords(kblock,1,1),coords(kblock,1,2),coords(kblock,1,3)/)
                    Xa = (/coords(kblock,2,1),coords(kblock,2,2),coords(kblock,2,3)/)
                    Xb = (/coords(kblock,3,1),coords(kblock,3,2),coords(kblock,3,3)/)
                    Xc = (/coords(kblock,4,1),coords(kblock,4,2),coords(kblock,4,3)/)
                elseif (i==2) then
                    Xp = (/coords(kblock,2,1),coords(kblock,2,2),coords(kblock,2,3)/)
                    Xa = (/coords(kblock,1,1),coords(kblock,1,2),coords(kblock,1,3)/)
                    Xb = (/coords(kblock,3,1),coords(kblock,3,2),coords(kblock,3,3)/)
                    Xc = (/coords(kblock,4,1),coords(kblock,4,2),coords(kblock,4,3)/)
                elseif (i==3) then
                    Xp = (/coords(kblock,3,1),coords(kblock,3,2),coords(kblock,3,3)/)
                    Xa = (/coords(kblock,1,1),coords(kblock,1,2),coords(kblock,1,3)/)
                    Xb = (/coords(kblock,2,1),coords(kblock,2,2),coords(kblock,2,3)/)

                    Xc = (/coords(kblock,4,1),coords(kblock,4,2),coords(kblock,4,3)/)
                elseif (i==4) then
                    Xp = (/coords(kblock,4,1),coords(kblock,4,2),coords(kblock,4,3)/)
                    Xa = (/coords(kblock,1,1),coords(kblock,1,2),coords(kblock,1,3)/)
                    Xb = (/coords(kblock,2,1),coords(kblock,2,2),coords(kblock,2,3)/)
                    Xc = (/coords(kblock,3,1),coords(kblock,3,2),coords(kblock,3,3)/)
                end if 
                Xba(1) = Xb(1)-Xa(1)
                Xba(2) = Xb(2)-Xa(2)
                Xba(3) = Xb(3)-Xa(3)
                Xca(1) = Xc(1)-Xa(1)
                Xca(2) = Xc(2)-Xa(2)
                Xca(3) = Xc(3)-Xa(3)
                Xpa(1) = Xp(1)-Xa(1)
                Xpa(2) = Xp(2)-Xa(2)
                Xpa(3) = Xp(3)-Xa(3)
                pNb = cross(Xba,Xca)
                pN = pNb/norm(pNb)
                pd = abs(dot(Xpa,pN))
                if (i==1) then
                    pd_min = pd
                else
                    if ( pd .lt. pd_min ) then
                        pd_min = pd
                    end if
                end if
                end do
                dtimeStable(kblock) = factorStable*(pd_min/cd)
                !dtimeStable(kblock) = 2.143102d-06
                
                !if (kblock==1) then
                !    write(*,*)"dtimeStable(kblock)"
                !    write(*,*)dtimeStable(kblock)
                !    write(*,*)"rhs(kblock=1,:)"
                !    write(*,*)rhs(kblock,:)
                !end if
                
            end do !----------------------------------------------------------------
            
            
            
        end if
            
    end if
                
    ! ============================================
    ! element outgoing note
    !write(*,*)" "
    !write(*,*)"UEL out - ELEMNO",ELEMNO
    !write(*,*)" "
    ! ============================================

return
end subroutine VUEL
!-------------------------------------------------------------------------------
!===============================================================================
