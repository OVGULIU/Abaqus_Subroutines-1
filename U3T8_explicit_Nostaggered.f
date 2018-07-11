        !DEC$ FREEFORM
        !===============================================================================
        ! DEVELOPER: EMMA GRIFFITHS & EDGAR HUSSER
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
            if (size(mat,1) == size(mat,2) .and. (size(mat,2) == 4) ) then
                    trace = mat(1,1)+ mat(2,2)+ mat(3,3)+mat(4,4)
            elseif (size(mat,1) == size(mat,2) .and. (size(mat,2) == 3) ) then
                    trace = mat(1,1)+ mat(2,2)+ mat(3,3)
            elseif (size(mat,1) == size(mat,2) .and. (size(mat,2) == 2) ) then
                    trace = mat(1,1)+ mat(2,2)
            else
                write(*,*)"Error in function `trace' - matrix is non-sqare!"
                CALL XPLB_EXIT
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
                write(*,*) "Error in function `det' - tensor is non-sqare!"
                CALL XPLB_EXIT
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
                write(*,*) "Error in function `inv' - tensor is non-sqare or larger then 3x3!"
                CALL XPLB_EXIT 
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
                do i= 1,size(vec1)
                    dot = dot + vec1(i)*vec2(i)
                end do
            else
                write(*,*) "Error in function `dot' - vectors have not the same length!"
                CALL XPLB_EXIT
            end if
    
        end function dot
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
                write(*,*) "error in `cross' - vector lengths are not the same or not given in 3D!"
                CALL XPLB_EXIT
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
                CALL XPLB_EXIT 
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
            write(*,*) "Error in function `dya' - vector lengths are not the same!"
                CALL XPLB_EXIT 
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
                write(*,*) "Dimension error in function `matvec' - dimension of vector must be consistent with the dimensions of the matrix"
                CALL XPLB_EXIT 
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
                write(*,*) "Dimension error in function `matmat' - matrix dimensions are not consistent or larger than 3x3"
                CALL XPLB_EXIT 
            end if

        end function matmat
        !===========================================================================
    
    end module functions
    !===============================================================================
    
!===============================================================================================================================================   
    module constitutive_relations

           implicit none

           !---------------------------------------------------------------------------
           public :: STRESSES_CONCEN_COUPLED,DIFFUSION_MECH_COUPLED
           
           contains
       
	!===========================================================================
            ! mechanical stresses resuling from coupled mechanical, diffusion and electric fields
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

!                S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID - (pEMCoup/pZ)*pQf*pID
		!write(*,*)"S",S
		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID + (1.d0/(pEPSILONZERO*pEPSILONR))*(dya(ElecDisp,ElecDisp) - 0.5d0*(dot(ElecDisp,ElecDisp))*pID)- (pEMCoup/pZ)*pQf*pID
!		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID
!		S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID + (1.d0/(pEPSILONZERO*pEPSILONR))*(dya(ElecDisp,ElecDisp) - 0.5d0*(dot(ElecDisp,ElecDisp))*pID)
            
                !---------------------------------------------------------------------------       
            
            end subroutine STRESSES_CONCEN_COUPLED

	!===========================================================================
            ! diffusive fluxes resuling from coupled mechanical, diffusion and electric fields
            subroutine DIFFUSION_MECH_COUPLED(J,pELECFIELD,Ee,pCo,gCo,pF,pZ,pEMCoup,pRTHETA,pDif,pStability)
       
               implicit none
                !---------------------------------------------------------------------------  
                ! declaration
         
                double precision, dimension(:), intent(inout) :: J
                double precision, dimension(:,:), intent(in) :: Ee
                double precision, dimension(:), intent(in) :: gCo, pELECFIELD, pStability
                double precision, intent(in) :: pCo,pF,pZ,pEMCoup,pRTHETA,pDif
                !---------------------------------------------------------------------------
		!write(*,*) "pCo", pCo
		!write(*,*)"pELECFIELD", pELECFIELD
		!write(*,*)"(pF*pZ)/(pRTHETA)*pDif*pCo*pELECFIELD", (pF*pZ)/(pRTHETA)*pDif*pCo*pELECFIELD
!		J = -pDif*gCo+(pF*pZ)/(pRTHETA)*pCo*pDif*pELECFIELD+(pF*pZ)/(pRTHETA)*pDif*pStability*gCo*pELECFIELD
!		J = -pDif*gCo+(pF*pZ)/(pRTHETA)*pDif*pELECFIELD
!		J = (pF*pZ)/(pRTHETA)*pDif*pELECFIELD
!		J = -pDif*gCo+pDif*pStability*gCo
!		J = -pDif*gCo-(pF*pZ)/(pRTHETA)*pDif*pStability*gCo*pELECFIELD
!		J = pDif*gCo
                !---------------------------------------------------------------------------       
            
            end subroutine DIFFUSION_MECH_COUPLED
    
        end module constitutive_relations

!===============================================================================================================================================
        !-------------------------------------------------------------------------------
        !
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

	kstep = i_Array(i_int_kStep)
	kInc = i_Array(i_int_kInc)
	
! ------ START OF THE ANALYSIS ------
	if (lOp .eq. j_int_StartAnalysis) then

! ------ Start of the step ------
	else if (lOp .eq. j_int_StartStep) then
! ------ Setup the increment ------
	else if (lOp .eq. j_int_SetupIncrement) then
! ------ Start of increment ------
	else if (lOp .eq. j_int_StartIncrement) then
		if (MOD(int(i_Array(i_int_kInc)),17000000).eq.0d0) then
			if (i_Array(i_int_kInc).ne.1) then
				write(*,*)
				write(*,*) "----------------- VEXTERNALDB at Increment:",int(i_Array(i_int_kInc)),"-----------------"
				write(*,*)
!				call EXECUTE_COMMAND_LINE('/home/cerecam/Desktop/GPG_Cube/PythonRunTest.sh', wait=.true.)
			end if
		end if
! ------ End of increment ------
	else if (lOp .eq. j_int_EndIncrement) then
! ------ End of the step ------
	else if (lOp .eq. j_int_EndStep) then
!			call EXECUTE_COMMAND_LINE('/home/cerecam/Desktop/VUEL_Simulations/VUEL_CUBE/Cube_Staggered_PNP_only/PythonRunTest.sh', wait=.true.)
! ------ End of the analysis ------
	else if (lOp .eq. j_int_EndAnalysis) then
!			write(*,*)
!			write(*,*) "----------------- VEXTERNALDB at Increment:",int(i_Array(i_int_kInc)),"-----------------",int(i_Array(i_int_kInc))
!			write(*,*)

	end if

	return
	end subroutine VEXTERNALDB

!===============================================================================================================================================
        !-------------------------------------------------------------------------------
        !
        ! USER SUBROUTINE - VUFIELD:UPDATES PREDEFINED FIELD
	SUBROUTINE VUFIELD(FIELD, NBLOCK, NFIELD, KFIELD, NCOMP, &
			   KSTEP, JFLAGS, JNODEID, TIME, &
			   COORDS, U, V, A)
	include 'vaba_param.inc'	
	! indices for the time array 
	integer, parameter :: i_ufld_Current = 1 
	integer, parameter :: i_ufld_Increment = 2 
	integer, parameter :: i_ufld_Period = 3 
	integer, parameter :: i_ufld_Total = 4

	! indices for the coordinate array COORDS 
	integer, parameter :: i_ufld_CoordX = 1 
	integer, parameter :: i_ufld_CoordY = 2 
	integer, parameter :: i_ufld_CoordZ = 3

	! indices for the displacement array U 
	integer, parameter :: i_ufld_SpaDisplX = 1 
	integer, parameter :: i_ufld_SpaDislY = 2 
	integer, parameter :: i_ufld_SpaDisplz = 3 
	integer, parameter :: i_ufld_RotDisplX = 4 
	integer, parameter :: i_ufld_RotDisplY = 5 
	integer, parameter :: i_ufld_RotDisplZ = 6 
	integer, parameter :: i_ufld_AcoPress = 7 
	integer, parameter :: i_ufld_Temp = 8

	!indices for velocity array V 
	integer, parameter :: i_ufld_SpaVelX = 1 
	integer, parameter :: i_ufld_SpaVelY = 2 
	integer, parameter :: i_ufld_SpaVelZ = 3 
	integer, parameter :: i_ufld_RotVelX = 4 
	integer, parameter :: i_ufld_RotVelY = 5 
	integer, parameter :: i_ufld_RotVelZ = 6 
	integer, parameter :: i_ufld_DAcoPress = 7 
	integer, parameter :: i_ufld_DTemp = 8

	! indicies for the acceleration array A 
	integer, parameter :: i_ufld_SpaAccelX = 1 
	integer, parameter :: i_ufld_SpaAccelY = 2 
	integer, parameter :: i_ufld_SpaAccelZ = 3 
	integer, parameter :: i_ufld_RotAccelX = 4 
	integer, parameter :: i_ufld_RotAccelY = 5 
	integer, parameter :: i_ufld_RotAccelZ = 6 
	integer, parameter :: i_ufld_DDAcoPress = 7 
	integer, parameter :: i_ufld_DDTemp = 8

	! indices for JFLAGS 
	integer, parameter :: i_ufld_kInc = 1 
	integer, parameter :: i_ufld_kPass = 2

	! Variables passed in for information
        integer, intent(in) :: NBLOCK
	integer, intent(in) :: NFIELD
	integer, intent(in) :: KFIELD
	integer, intent(in) :: NCOMP
	integer, intent(in) :: KSTEP
	integer, intent(in), dimension(i_ufld_kPass) :: JFLAGS
	integer, intent(in), dimension(NBLOCK) :: JNODEID
	double precision, intent(in), dimension(4) :: TIME
	double precision, intent(in), dimension(3,NBLOCK) :: COORDS
	double precision, intent(in), dimension(8,NBLOCK) :: U,V,A
	double precision, dimension(NBLOCK) :: data_arr
	

	! Dimensioned arrays
	double precision, intent(inout), dimension(NBLOCK,NCOMP,NFIELD) :: FIELD
!	write(*,*) "TIME(i_ufld_Current)",TIME(i_ufld_Current)
!	write(*,*) "MOD(TIME(i_ufld_Current),0.05d0)",MOD(TIME(i_ufld_Current),0.05d0)
!	write(*,*) "MOD(TIME(i_ufld_Current),0.05d0).eq.0.0d0", MOD(TIME(i_ufld_Current),0.05d0).eq.0.0d0
!	write(*,*) "TIME(i_ufld_Current).ne.0",TIME(i_ufld_Current).ne.0
!	write(*,*)

!			write(*,*) "NCOMP", NCOMP, "NFIELD", NFIELD
!			write(*,*) "NBLOCK", NBLOCK
!			write(*,*) "KSTEP", KSTEP
!			write(*,*) "JFLAGS", JFLAGS
!			write(*,*) "JNODEID(1)", JNODEID(1)
!			write(*,*) "JNODEID(15)", JNODEID(15)
!			write(*,*) "JNODEID(55)", JNODEID(55)
!			write(*,*) "JFLAGS",JFLAGS
!			write(*,*) "U(all,1)", U(:,kblock)
!			write(*,*) "FIELD",FIELD(:,NCOMP,NFIELD)
!	if (MOD(JFLAGS(i_ufld_kInc),1).eq.0d0) then
			if (JFLAGS(i_ufld_kInc).eq.1d0 .OR. JFLAGS(i_ufld_kInc).eq.0d0) then				
			elseif (MOD(int(JFLAGS(i_ufld_kInc)),17000000).eq.0d0) then
				write(*,*) "----------------- VUFIELD at increment:",JFLAGS(i_ufld_kInc)," -----------------"	
				open(unit=105, file='/home/cerecam/Desktop/GPG_Cube/ElecPotentialsSandwich.csv',status='old')!
				READ(105,*) data_arr
				do kblock=1,nblock
					FIELD(kblock,NCOMP,NFIELD) = data_arr(kblock)
				end do
!				write(*,*) "********* New field written **********"
!				write(*,*)
				close(105)
			end if
!	end if
!		write(*,*) "VUFIELD"

	return
	end subroutine VUFIELD

!===============================================================================================================================================
        !
        !-------------------------------------------------------------------------------
        !
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
            
            implicit none ! use double=both and output_precision=full
            !include 'vaba_param.inc'
            ! if implicit none is not used and 'vaba_param.inc' is included - implicit declaration of variables:
            ! all variables for letters a-h, o-z is real and i-n are integers (i,j,k,l,m,n).
            !-------------------------------------------------------------------------------   
            ! declaration
            
            ! parameters - numbers
            double precision, parameter :: zero=0.d0
            double precision, parameter :: half=0.5d0
            double precision, parameter :: one=1.d0
            double precision, parameter :: two=2.d0
            double precision, parameter :: three=3.d0
            double precision, parameter :: four=4.d0
            double precision, parameter :: six=6.d0
            double precision, parameter :: eight=8.d0
            double precision, parameter :: Abcissa=SQRT(1.d0/3.d0)
            double precision, parameter :: factorStable=0.99d0
            double precision, parameter :: pi=3.1415926535897932
            
            ! parameters - problem specification
            integer, parameter :: iGP=8		!Number of Gauss Points
            integer, parameter :: iCORD=3	!Degrees of freedom (mechanical)
            integer, parameter :: ICORDTOTAL=4  !Degrees of freedom (total per node (3 mech; 1 temp))
            integer, parameter :: iNODE=8	!Number of nodes per element

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
            
            !  predefined parameters - procedure flags (Explicit Dynamic)
            integer, parameter :: jDynExplicit = 17
	    !  predefined parameters - procedure flags
            integer, parameter :: jDynThermalStressExplicit = 74
            
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
            integer, intent(in ) :: kstep 	                        ! current step number
            integer, intent(in ) :: kinc 	                        ! current increment number
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
            double precision, intent(inout), dimension(nblock,ndofel,ndofel) :: amass  	! element contribution to mass matrix (symmetric)
            double precision, intent(inout), dimension(nblock) :: dtimeStable  	    	! scalar value, upper limit of time increment (Courant condition)
            double precision, intent(inout), dimension(nblock,nsvars) :: svars  	! solution-depending state variables
            double precision, intent(inout), dimension(nblock,nElEnergy) :: energy	! element energy quantities 
            
            ! user variables
            ! ==================
            
            double precision :: pGPCORD(iGP,iCORD)  ! integration point coordinates
            double precision :: pQUAD               ! factor for quadrature rule
            double precision :: pWT(iGP)            ! integration point weights
            double precision :: pNN(iGP,iNODE)      ! shape function matrix (interpolation matrix)
            double precision :: pID(iCORD,iCORD)    ! identity tensor
            
            double precision :: amass_row_sum       ! sum of array row
            double precision :: cd,pd,pd_min,cdT, pd_min_Jelem
            double precision ::	Mechtime,Thermaltime,TimeMin
            ! include others

            !===============================================================================

            double precision :: xi1,xi2,xi3         ! natural coordinates
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
            double precision :: Uarray(iCORD)       ! nodal displacement in each dimension
            double precision :: CoNODE(iNODE)       ! nodal concentration 1.4498201477996281E-007s
            double precision :: Co           ! Concentration inetrpolated in element
            double precision :: pCo          ! Concentration mobile interpolated in element
            double precision :: pQf       ! Free charge density
            double precision :: phi          ! nodal Electrical potentials
            double precision :: gCo(iCORD)          ! Grad of concentrations
            double precision :: gCo2(iCORD)          ! Grad of concentrations
            double precision :: H(iCORD,iCORD)      ! displacement gradient
            double precision :: F(iCORD,iCORD)      ! deformation gradient
            double precision :: Ee(iCORD,iCORD)     ! small /Gree-Lagrange strain tensor
            double precision :: SIG(iCORD,iCORD)    ! Cauchy stress tensor
            double precision :: S(iCORD,iCORD)      ! 2. PK
            double precision :: P(iCORD,iCORD)      ! 1. PK
            double precision :: SIGt(iCORD)    ! Cauchy stress tensor
            double precision :: St(iCORD)      ! 2. PK
            double precision :: Pt(iCORD)      ! 1. PK
            double precision :: pELECFIELD(iCORD)    ! Electric Field vector from grad(phi)
            double precision :: ELECDISP(iCORD)    ! Electric Displacement vector calculated from elec field vector
            double precision :: J(iCORD)    ! flux vector for diffusion problem
            double precision :: VonMisS    ! Von Mises stress
            double precision :: VonMisE    ! Von Mises Strain
            integer :: dofni(iCORD)        ! current dof
            integer :: dofniT        	   ! Thermal current dof
            integer :: dofnj(iCORD)        ! current dof
            integer :: dofnjT        	   ! Thermal current dof

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
            double precision :: Csat		! Saturated concentration
            double precision :: pVolume		! Volume of element
            double precision :: pSED		! SED
            double precision :: pa1		! stability term RHS
            double precision :: pa2		! stability term AMASS
            double precision :: Pe		! Peclet number
            double precision :: Courant		! Courant number
            double precision :: Sigma_K     ! Stability parameter [ref. Tu et al, Compt Physics Comm]
            double precision :: pA_VECTOR(iCORD)	! Stabilization vector [ref. Tu et al, Compt Physics Comm]
            double precision :: V_DASH		! Additive diffusion
            double precision :: elesize_sum	! Element size intermediate term
            double precision :: Elesize		! Element size in veocity field direction (Radius of equivalent volume circle)
            double precision :: pCoTotal	! integration within domain
            double precision :: pStability(iCORD)	! stability
                
            
            ! integer
            integer :: kblock,ip,nn,ni,nj,i


            integer :: CordRange(iCORD) = (/(i, i=1,iCORD, 1)/)
            integer :: NODERange(iNODE) = (/(i, i=1,iNODE, 1)/)
            
            !-------------------------------------------------------------------------------       

!===============================================================================================================================================

                      
        if (jtype.eq.38) then        
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
            pa2 = props(13)
            
            pGM  = half*pEM/(one+pNU)
            pLAM = (pEM*pNU)/((one+pNU)*(one-two*pNU))    
            	
            ! integration point coordinates and weights --------------------------------
            if (iGP==8) then ! HUGHES - The Finite Element Method 1987 (p. 174)
    
                pGPCORD(1,:) = (/ -Abcissa, -Abcissa, -Abcissa /)
                pGPCORD(2,:) = (/  Abcissa, -Abcissa, -Abcissa /)
                pGPCORD(3,:) = (/  Abcissa,  Abcissa, -Abcissa /)
                pGPCORD(4,:) = (/ -Abcissa,  Abcissa, -Abcissa /)
                pGPCORD(5,:) = (/ -Abcissa, -Abcissa,  Abcissa /)
                pGPCORD(6,:) = (/  Abcissa, -Abcissa,  Abcissa /)
                pGPCORD(7,:) = (/  Abcissa,  Abcissa,  Abcissa /)
                pGPCORD(8,:) = (/ -Abcissa,  Abcissa,  Abcissa /)
                    
                pWT(:) = (/ one, one, one, one, one, one, one, one /)
                    
                pQUAD = one
                    
            else
    
                write(*,*) "Error in computation of integration point coordinates. The number of IPs does not conform with the element type (4 node tetrahedral element with 4 integratin points)."
                CALL XPLB_EXIT
            end if
                                
            ! shape function for each ip -----------------------------------------------
            do ip=1,iGP
    
                xi1=pGPCORD(ip,1)
                xi2=pGPCORD(ip,2)
                xi3=pGPCORD(ip,3)
    
                if (iNODE==8) then ! cf. WRIGGERS - Nonlinear Finite Element3 Methods 2008 (p. 120)
    
                    pNN(ip,1) = one/eight*(1-xi1)*(1-xi2)*(1-xi3)
                    pNN(ip,2) = one/eight*(1+xi1)*(1-xi2)*(1-xi3)
                    pNN(ip,3) = one/eight*(1+xi1)*(1+xi2)*(1-xi3)
                    pNN(ip,4) = one/eight*(1-xi1)*(1+xi2)*(1-xi3)
                    pNN(ip,5) = one/eight*(1-xi1)*(1-xi2)*(1+xi3)
                    pNN(ip,6) = one/eight*(1+xi1)*(1-xi2)*(1+xi3)
                    pNN(ip,7) = one/eight*(1+xi1)*(1+xi2)*(1+xi3)
                    pNN(ip,8) = one/eight*(1-xi1)*(1+xi2)*(1+xi3)
    
                else
                    write(*,*) "Error in computation of shape functions. The number of nodes does not conform with the element type (4 node tetrahedral element)."
                    CALL XPLB_EXIT
                end if
    
            end do
            !---------------------------------------------------------------------------
                
            ! identity matrix ----------------------------------------------------------
            pID=zero
            forall(i=1:iCORD) pID(i,i)=one
            !---------------------------------------------------------------------------

    !===============================================================================================================================================
            
            ! loop over element block
            do kblock=1,nblock ! ---------------------------------------------------
            
                ! loop over all integration points (computation of FE variables)
                do ip=1,iGP ! ------------------------------------------------------
        
        
                    ! natural coordinates of current ip
                    xi1 = pGPCORD(ip,1)
                    xi2 = pGPCORD(ip,2)
                    xi3 = pGPCORD(ip,3)
        
                    ! coordinate vectors of current element
                    X1 = coords(kblock,:,1)
                    X2 = coords(kblock,:,2)
                    X3 = coords(kblock,:,3)
                    
                    ! --------------------------------------------------------------
                    if (iNODE==8) then
                        
                        ! derivatives of shape functions with respect to natural coordinates                
                        dNdXi1(1) = -one/eight*(1-xi2)*(1-xi3)
                        dNdXi2(1) = -one/eight*(1-xi1)*(1-xi3)
                        dNdXi3(1) = -one/eight*(1-xi1)*(1-xi2)
        
                        dNdXi1(2) =  one/eight*(1-xi2)*(1-xi3)
                        dNdXi2(2) =  -one/eight*(1+xi1)*(1-xi3)
                        dNdXi3(2) =  -one/eight*(1+xi1)*(1-xi2)
        
                        dNdXi1(3) =  one/eight*(1+xi2)*(1-xi3)
                        dNdXi2(3) =  one/eight*(1+xi1)*(1-xi3)
                        dNdXi3(3) =  -one/eight*(1+xi1)*(1+xi2)
        
                        dNdXi1(4) =  -one/eight*(1+xi2)*(1-xi3)
                        dNdXi2(4) =  one/eight*(1-xi1)*(1-xi3)
                        dNdXi3(4) =  -one/eight*(1-xi1)*(1+xi2)
        
                        dNdXi1(5) =  -one/eight*(1-xi2)*(1+xi3)
                        dNdXi2(5) =  -one/eight*(1-xi1)*(1+xi3)
                        dNdXi3(5) =  one/eight*(1-xi1)*(1-xi2)
        
                        dNdXi1(6) =  one/eight*(1-xi2)*(1+xi3)
                        dNdXi2(6) =  -one/eight*(1+xi1)*(1+xi3)
                        dNdXi3(6) =  one/eight*(1+xi1)*(1-xi2)
        
                        dNdXi1(7) =  one/eight*(1+xi2)*(1+xi3)
                        dNdXi2(7) =  one/eight*(1+xi1)*(1+xi3)
                        dNdXi3(7) =  one/eight*(1+xi1)*(1+xi2)
        
                        dNdXi1(8) =  -one/eight*(1+xi2)*(1+xi3)
                        dNdXi2(8) =  one/eight*(1-xi1)*(1+xi3)
                        dNdXi3(8) =  one/eight*(1-xi1)*(1+xi2)
                     
                    else  
                        write(*,*) "Error in computation of shape function derivatives. The number of nodes does not conform with the element type (4 node tetrahedral element)." 
                        CALL XPLB_EXIT    
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
                    JJ(1,:) = (/dX1dxi1, dX2dxi1, dX3dxi1/)
                    JJ(2,:) = (/dX1dxi2, dX2dxi2, dX3dxi2/)
                    JJ(3,:) = (/dX1dxi3, dX2dxi3, dX3dxi3/)
                    
                    detJ(ip) = det(JJ)
                    
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
    
                    end do !----------------nn-loop --------------------
    
               
                end do ! -------------------ip-loop------------------------------- 
!    !===================================================================================================================================
!    !----------------------------------------------ELEMENT LENGTH CALCULATION----------------------------------------------------
!    !===================================================================================================================================
                
                pVolume = detJ(1)+detJ(2)+detJ(3)+detJ(4)+detJ(5)+detJ(6)+detJ(7)+detJ(8)
                pd_min = pvolume**(one/three)
                cd = sqrt( (pEM*(one-pNU))/(pRHO*(one+pNU)*(one-two*pNU)) )
                cdT = (pDif)
                Mechtime = (pd_min/cd)
                Thermaltime = (pd_min*pd_min)/(2*cdT)
                TimeMin = minval( (/Mechtime,Thermaltime/) )
                dtimeStable(kblock) = factorStable*TimeMin		
!    !===================================================================================================================================
!    !--------------------------------------------------------RHS CALCULATION------------------------------------------------------------
!    !===================================================================================================================================
    
                energy(kblock,iElIe)=zero
                energy(kblock,iElTh)=zero
<<<<<<< HEAD
                rhs(kblock,1:ndofel)=zero
                !energy(kblock,iElKe)=zero
                    
                if (lflags(iOpCode).eq.jIntForceAndDtStable) then
=======
                energy(kblock,iElKe)=zero
                    
!                if (lflags(iOpCode).eq.jIntForceAndDtStable) then
>>>>>>> 6b7c65d97ff645583f3fafcbd374ed3881e23f10
                    do ip=1,iGP ! ---------------------- loop over all integration points (computation of residuum)--------------------------------
                    
                        H = zero
                        pQf = zero
                        pCo = zero
                        gCo = zero
                        gCo2 = zero
                        pELECFIELD= zero
                        Elesize = zero
                        do ni=1,iNODE
                            dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
                            Uarray(1:iCORD)=u(kblock, dofni)
        
                            ! displacement gradient
                            H = H + dya(Uarray,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
            
                            ! Concentration and Concentration gradient
                            CoNODE(ni) = u(kblock, (iCORDTOTAL*ni))
                            gCo = gCo + ( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )*CoNODE(ni)
                            
                            !Electric Field
                            pELECFIELD = pELECFIELD - ((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))*predef(kblock,ni,2,1)	
                        end do
<<<<<<< HEAD
=======
                        
                        Elesize = pd_min
                        Pe = NORM((pZ*pF/pRTHETA*pELECFIELD*(Elesize)))/2
                        Courant = NORM((pDif*pZ*pF/pRTHETA*pELECFIELD))*dtimeCur/Elesize
        !                if (Pe>1.0) then
        !                end if
        !                if (Courant>0.1)   then
        !                    write(*,*) "Courant number: ", Courant, "at ", jElem(kblock)
        !                end if
                        
                        pCo = dot(pNN(ip,:),CoNODE)
>>>>>>> 6b7c65d97ff645583f3fafcbd374ed3881e23f10
               
                        ! small strain tensor
                        Ee = half*(transpose(H) + H)
        
<<<<<<< HEAD
                        ! Elemental concentration and charge carrier density
                        pCo = dot(pNN(ip,:),CoNODE)
                        pQf = pF*((pZ*pCo)+(cSat*(1.d0)))
        
                        ! Electrical Displacement given by -(minus).epsilon0.epsilonr.Elecfield
                        ElecDisp = pEPSILONZERO*pEPSILONR*pELECFIELD
                        
                        ! Internal energy Calculation							
                        pSED = ( ( pGM*ddot(Ee,Ee)+half*pLAM*trace(Ee)*trace(Ee)) )            
                        energy(kblock,iElIe)= energy(kblock,iElIe) + (detJ(ip))*pSED
            
                        ! Stability calculations
                        Elesize = pd_min
                        Pe = NORM((pZ*pF/pRTHETA*pELECFIELD*(Elesize)))/2
                        Courant = NORM((pDif*pZ*pF/pRTHETA*pELECFIELD))*dtimeCur/Elesize
                        pA_Vector = pDif*pZ*pF/pRTHETA*pELECFIELD                    
                        sigma_k = (Elesize/(2*NORM(pA_Vector)))*Pe
                        
!                        if (Pe>1.0) then
!                        end if
!                        if (Courant>0.1)   then
!                            write(*,*) "Courant number: ", Courant, "at ", jElem(kblock)
!                        end if
        
                        do ni=1,iNODE !-----------------------------loop-i--------------
=======
                        ! Electrical Displacement given by -(minus)X epsilon0 Xepsilonr XElecfield
                        ElecDisp = pEPSILONZERO*pEPSILONR*pELECFIELD
        
                        rhs(kblock,1:ndofel)=zero
    !										
                        pSED = ( ( pGM*ddot(Ee,Ee)+half*pLAM*trace(Ee)*trace(Ee)) )
            
                        energy(kblock,iElIe)= energy(kblock,iElIe) + (detJ(ip))*pSED
            
                        pQf = pF*((pZ*pCo)+(cSat*(1.d0)))
                        
                        pA_Vector = pDif*pZ*pF/pRTHETA*pELECFIELD                    
                        sigma_k = (Elesize/(2*NORM(pA_Vector)))*Pe
                        do ni=1,iNODE !-----------------------------loop-i--------------
                            pQf = pF*((pZ*pCo)+(cSat*(1.d0)))
>>>>>>> 6b7c65d97ff645583f3fafcbd374ed3881e23f10
    
                            call STRESSES_CONCEN_COUPLED(S,Ee,ElecDisp,pQf,pID,pGM,pLAM,pEPSILONZERO,pEPSILONR,pEMCoup, pZ, cSat)
                            
                            dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
                            dofniT = iCORDTOTAL*ni
                        
                    !--------------------------------------Displacement RHS--------------------------------------
                            rhs(kblock,dofni) = rhs(kblock,dofni) 	+ pQUAD*pWT(ip)*detJ(ip)*(matvec(S,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))) 
    !											- pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),((pEMCoup/pZ)*pQf*pID))
                    !--------------------------------------Concentration RHS--------------------------------------
                            rhs(kblock,dofniT) = rhs(kblock,dofniT) &
                        - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(-gCo)) &
                        - pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(((pF*pZ)/(pRTHETA)*pNN(ip,ni)*pCo*pELECFIELD))) &
<<<<<<< HEAD
                        + (pF*pZ)/(pRTHETA)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,(sigma_k*pA_Vector)*pa1)))
=======
                        - (pF*pZ)/(pRTHETA)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,(sigma_k*pA_Vector)*pa1)))
>>>>>>> 6b7c65d97ff645583f3fafcbd374ed3881e23f10
    !					+ pDif*(pF*pZ)/(pRTHETA)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,((/1.0d0, 1.0d0, 1.0d0/)*pa1))))
    !                   
    !         				rhs(kblock,dofniT) = rhs(kblock,dofniT) &
    !					- pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(-pDif*gCo))
        
    !         				rhs(kblock,dofniT) = rhs(kblock,dofniT) &
    !					- pQUAD*pWT(ip)*detJ(ip)*dot( (/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),((pDif*(pF*pZ)/(pRTHETA)*pNN(ip,ni)*pCo*pELECFIELD)))&
    !					+ pDif*(pF*pZ)/(pRTHETA)*pQUAD*pWT(ip)*detJ(ip)*dot((/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/),(gCo*dot(pELECFIELD,((/1.0d0, 1.0d0, 1.0d0/)*pa))))
    
                    !--------------------------------------Thermal Energy--------------------------------------
                            energy(kblock,iElTh)= energy(kblock,iElTh) + (pDif*pNN(ip,ni)*u(kblock,dofniT))
                            
                    !--------------------------------------Kinetic Energy--------------------------------------
                            !do nj=1,iNODE !-----------------------------loop-i--------------
                            !	dofnj(1:iCORD) = 1+iCORDTOTAL*((nj*1)-1)+(CordRange-1)
                            !	energy(kblock,iElKe)= energy(kblock,iElKe) + half*dot(v(kblock,dofni),matvec(amass(kblock,dofni,dofnj),v(kblock,dofnj)))
                            !end do !------------------------------end-loop-nj----------------
                        end do !------------------------------end-loop-ni----------------
    
                    end do ! -------------------- ip loop  ------------------------
                else if (lflags(iOpCode).eq.jMassCalc) then
                    !amass(kblock,:,:)= zero
                    do ip=1,iGP ! ---------------------- loop over all integration points (computation of mass matrix)--------------------------------
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
                                                            + (1/pDif)*pQUAD*pWT(ip)*detJ(ip)*pNN(ip,ni)*pNN(ip,nj)
                                        
                            end do !--------------------------end-loop-j----------------
                        
                        end do !------------------------------end-loop-i----------------
                        ! mass lumping
                        
                    end do ! --------------------- end-loop-ip -----------------------
                    
                    do i=1,ndofel
                        amass_row_sum = sum(amass(kblock,i,:))
                        amass(kblock,i,:) = zero
                        amass(kblock,i,i) = amass_row_sum
                    end do
                end if ! ----------------------------jMassCalc--------------------
    
                               
    !===================================================================================================================================
    !----------------------------------------------STABLE TIME INCREMENT CALCULATION----------------------------------------------------
    !===================================================================================================================================
    
                
                            
                           
    !===============================================================================================================================================
                                                    
            end do !-----------------------nblock-----------------------------------
        
        end if ! ------------------------ type.eq.48 loop ---------------------------
    
    return
    end subroutine VUEL
        !-------------------------------------------------------------------------------
        !===============================================================================
    			
