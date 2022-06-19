!
!    Copyright (C) 2017, P. R. Wiecha, A. Arbouet, C. Girard
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!    
!
!*********************************************************************
! Original Author: Christian Girard 
! CEMES/CNRS 
!
! Date: 24 January 1996 
! Modified: April 1997 
! Including the renormalisation terms
!
! ***********
! Modified 2010-2012 by A. Arbouet, R. Marty
! CEMES/CNRS
!   - openMP parallelization of matrix setup, Dyson sequence
!
! ***********
! Modified 2014-2016 by P.R. Wiecha
! CEMES/CNRS
!   - separated matrix setup / Dyson sequence
!
!*********************************************************************




!*********************************************************************
!--------------------------------------------------------------------*
!         Implements the the side-by-side matrix "S" containing      *
!         the coupled free-space and surface propagators             *
!         for the system (XM,YM,ZM).                                 *
!         The inverse of this matrix is the generalized propagator K *
!         with the property: K*E0 = E                                *
!--------------------------------------------------------------------*
!*********************************************************************
      SUBROUTINE SETUPMATRIX(ALAMBDA,SPACE, &
	          & XM,YM,ZM,  &
	          & CAP, &
	          & cn1,cn2,cn3, &
	          & CNORM, TIMESALPHA, PROPAGATOR, &
	          & NMAX, &
	          & CAM0)
! 	          NTHREADS, &
      
      !f2py threadsafe
      !$ use OMP_LIB
      USE OMP_LIB
      USE PYGDMPRECISION, ONLY : dp
      implicit none
      
      
!**********************************************************************
!*      CONFIGURATION
!**********************************************************************
      ! Fixed allocations
      real(dp), parameter :: Pi=3.141592654_dp
      complex(dp), parameter :: C0=(0._dp,0._dp)
      complex(dp), parameter :: CUN=(1._dp,0._dp)
      complex(dp), parameter :: CIM=(0._dp,1._dp)
     
     
      integer :: I,J,II,JJ
      integer, intent(in) :: NMAX 
      integer, intent(in) :: TIMESALPHA
!       integer, intent(in) :: NTHREADS
      complex(dp),   intent(in) :: cn1,cn2,cn3
      real(dp), intent(in) :: ALAMBDA
      real(dp), intent(in) :: SPACE
      real(dp), intent(in) :: XM(NMAX),YM(NMAX),ZM(NMAX)
      complex(dp), intent(in) :: CAP(NMAX)  ! complex polarizabilities of each meshpoint
      complex(dp), intent(in) :: CNORM  ! mesh-type dependent normalization factor
      
      !*** --- PROPAGATOR TO USE: --- ***
      !         - 1: 'simple 123'
      !         - 2: 'full 123'
      integer, intent(in) :: PROPAGATOR
      !**********************************
      
      complex(dp), intent(out) :: CAM0(3*NMAX,3*NMAX)
     
     
     
!**********************************************************************
!*      DECLARATION
!**********************************************************************
!       integer :: NTHREADSUSE
      real(dp) :: AK00
      real(dp) :: XI,YI,ZI,XJ,YJ,ZJ,ZS,XD,YD,ZD
      complex(dp) :: CHXX,CHXY,CHXZ,CHYX,CHYY,CHYZ,CHZX,CHZY,CHZZ
      complex(dp) :: SXX,SXY,SXZ,SYY,SYZ,SZZ
      complex(dp) :: CTXX,CTXY,CTXZ,CTYY,CTYZ,CTZZ
      complex(dp) :: CSXX,CSXY,CSXZ,CSYY,CSYZ,CSZZ!,CSYX,CSZX,CSZY
      real(dp) :: TXX,TXY,TXZ,TYY,TYZ,TZZ
      real(dp) :: TXX2,TXY2,TXZ2,TYY2,TYZ2,TZZ2
      
      complex(dp) :: CAXX(NMAX), CAYY(NMAX), CAZY(NMAX)
      complex(dp) :: CAZZ(NMAX), CAXY(NMAX), CAXZ(NMAX)
      complex(dp) :: CAYX(NMAX), CAYZ(NMAX), CAZX(NMAX)
     
      
      
      
      
!       IF (NTHREADS==-1) THEN
!         NTHREADSUSE = OMP_GET_MAX_THREADS()
!       ELSE
!         NTHREADSUSE = NTHREADS
!       ENDIF
!       CALL omp_set_num_threads(NTHREADSUSE)
 
!*********************************************************************
!*          INITIALIZATION
!*********************************************************************
      AK00=2._dp*Pi/ALAMBDA
     
      if (TIMESALPHA==1) then
!-----(1) Definition of the anisotropic polarizability of each -------
!---------------------discretized  cell ------------------------------
        DO J=1,NMAX
            CAXX(J)=CAP(J)
            CAYY(J)=CAP(J)
            CAZZ(J)=CAP(J)
            CAXY(J)=C0
            CAYX(J)=C0
            CAXZ(J)=C0
            CAZX(J)=C0
            CAZY(J)=C0
            CAYZ(J)=C0
        ENDDO
     else
     !! if demanded, don't multiply by polarizability (for python-treated inversion)
     !! --> use unitary tensor for polarizability
        DO J=1,NMAX
            CAXX(J)=CUN
            CAYY(J)=CUN
            CAZZ(J)=CUN
            CAXY(J)=C0
            CAYX(J)=C0
            CAXZ(J)=C0
            CAZX(J)=C0
            CAZY(J)=C0
            CAYZ(J)=C0
        ENDDO
     endif


!---------Initialization of the matrix composed with ----------------
!-------------------site-to-site propagators-------------------------
!--------------------------AM0(I,J)----------------------------------
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(XJ,YJ,ZJ,JJ, XI,YI,ZI,II, XD,YD,ZD,ZS, &
!$OMP&              SXX,SYY,SZZ,SXY,SXZ,SYZ,  &
!$OMP&              CHXX,CHYY,CHZZ,CHXY,CHYX,CHYZ,CHZY,CHXZ,CHZX, &
!$OMP&              TXX,TYY,TZZ,TXY,TXZ,TYZ,TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2, &
!$OMP&              CTXX,CTYY,CTZZ,CTXY,CTXZ,CTYZ, &
!$OMP&              CSXX,CSYY,CSZZ,CSXY,CSXZ,CSYZ &
!$OMP&              )

!$OMP DO
      DO J=1,NMAX
        
        XJ=XM(J)
        YJ=YM(J)
        ZJ=ZM(J)
        JJ=3*(J-1)

        DO I=1,NMAX
            XI=XM(I)
            YI=YM(I)
            ZI=ZM(I)
            II=3*(I-1)

!C*************************************************************************
!C****** simple model: vacuum propagator + 1/2/3-surface, material in center medium
!C*************************************************************************
            if (PROPAGATOR==1) then
            ! --- "self"-term
            IF(I.EQ.J)THEN
                XD=0._dp
                YD=0._dp
                ZS=ZI+ZJ
                
                !! --- Surface propagator
                CALL PROPAS(XD,YD,ZS,SXX,SYY,SZZ,SXY,SXZ,SYZ,SPACE, &
                          & cn1,cn2,cn3) 
                
                !! Setup matrix G. Final matrix to inverse will be: (1 - G.alpha)
                CHXX = (SXX+CNORM)
                CHYY = (SYY+CNORM)
                CHZZ = (SZZ+CNORM)
                     
                CHXY =   SXY*CUN
                CHYX =   SXY*CUN
                CHYZ =   SYZ*CUN
                CHZY = (-1._dp*SYZ*CUN)
                CHXZ =   SXZ*CUN
                CHZX = (-1._dp*SXZ*CUN)
                
            ! --- Couple 2 different dipoles
            ELSE
                    XD=XI-XJ
                    YD=YI-YJ
                    ZD=ZI-ZJ
                    ZS=ZI+ZJ
                
!                     write (*,*) "prop1"
                    !! --- Vacuum propagator
                    CALL PROPA0(AK00,XD,YD,ZD,TXX,TYY,TZZ,&
                            & TXY,TXZ,TYZ,TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2,&
                            & cn2)
                    
                    CTXX=TXX*CUN+TXX2*CIM
                    CTYY=TYY*CUN+TYY2*CIM
                    CTZZ=TZZ*CUN+TZZ2*CIM
                    CTXY=TXY*CUN+TXY2*CIM
                    CTXZ=TXZ*CUN+TXZ2*CIM
                    CTYZ=TYZ*CUN+TYZ2*CIM
                    
                    !! --- Surface propagator
                    CALL PROPAS(XD,YD,ZS,SXX,SYY,SZZ,SXY,SXZ,SYZ,SPACE,&
                            & cn1,cn2,cn3)
                    
                    CSXX=SXX*CUN
                    CSYY=SYY*CUN
                    CSZZ=SZZ*CUN
                    CSXY=SXY*CUN
                    CSXZ=SXZ*CUN
                    CSYZ=SYZ*CUN
                !------------------------------------------------------------------
                    CHXX = (CTXX+CSXX)
                    CHYY = (CTYY+CSYY)
                    CHZZ = (CTZZ+CSZZ)
                    CHXY = (CTXY+CSXY)
                    CHXZ = (CTXZ+CSXZ)
                    CHYZ = (CTYZ+CSYZ)
                    CHYX = (CTXY+CSXY)
                    CHZY = (CTYZ-CSYZ)
                    CHZX = (CTXZ-CSXZ)
                endif
                
! !C*************************************************************************
! !C****** magnetic propagator (in vacuum)
! !C*************************************************************************
                elseif (PROPAGATOR==2) then
                ! --- self-term
                IF(I.EQ.J)THEN
                ! --- "self"-term
                CHXX = CNORM
                CHYY = CNORM
                CHZZ = CNORM
                     
                CHXY = C0
                CHYX = C0
                CHYZ = C0
                CHZY = C0
                CHXZ = C0
                CHZX = C0
                
                ELSE
                ! --- Couple 2 different dipoles
!                   write (*,*) "prop2"
                    CALL ELECMAG0_FULL(AK00, &
                                        & XJ,YJ,ZJ, XI,YI,ZI, &
                                        & CHXX, CHXY, CHXZ, &
                                        & CHYX, CHYY, CHYZ, &
                                        & CHZX, CHZY, CHZZ, &
                                        & cn1,cn2,cn3)
                endif
            ENDIF
            
            
            ! multiply by anisotropic polarizability tensor
            ! "-" for  (1 - G*alpha); below treat only diagonal elements
            CAM0(II+1,JJ+1) = (CHXX*CAXX(J) + CHXY*CAYX(J) + CHXZ*CAZX(J))
            CAM0(II+2,JJ+1) = - (CHYX*CAXX(J) + CHYY*CAYX(J) + CHYZ*CAZX(J))
            CAM0(II+3,JJ+1) = - (CHZX*CAXX(J) + CHZY*CAYX(J) + CHZZ*CAZX(J))
            
            CAM0(II+1,JJ+2) = - (CHXX*CAXY(J) + CHXY*CAYY(J) + CHXZ*CAZY(J))
            CAM0(II+2,JJ+2) = (CHYX*CAXY(J) + CHYY*CAYY(J) + CHYZ*CAZY(J))
            CAM0(II+3,JJ+2) = - (CHZX*CAXY(J) + CHZY*CAYY(J) + CHZZ*CAZY(J))
            
            CAM0(II+1,JJ+3) = - (CHXX*CAXZ(J) + CHXY*CAYZ(J) + CHXZ*CAZZ(J))
            CAM0(II+2,JJ+3) = - (CHYX*CAXZ(J) + CHYY*CAYZ(J) + CHYZ*CAZZ(J))
            CAM0(II+3,JJ+3) = (CHZX*CAXZ(J) + CHZY*CAYZ(J) + CHZZ*CAZZ(J))
            
            
            ! Matrix to inverse: CAM0 = (1-G.alpha) 
            ! -->Substract from unitary tensor
            IF(I.EQ.J)THEN
                CAM0(II+1,JJ+1) = 1_dp - CAM0(II+1,JJ+1)
                CAM0(II+2,JJ+2) = 1_dp - CAM0(II+2,JJ+2)
                CAM0(II+3,JJ+3) = 1_dp - CAM0(II+3,JJ+3)
            ELSE
                CAM0(II+1,JJ+1) = - CAM0(II+1,JJ+1)
                CAM0(II+2,JJ+2) = - CAM0(II+2,JJ+2)
                CAM0(II+3,JJ+3) = - CAM0(II+3,JJ+3)
            ENDIF
            
            
! ------------ sub-tensor components in CAM0: ------------
!             CAM0(II+1,JJ+1) --> XX
!             CAM0(II+2,JJ+1) --> YX
!             CAM0(II+3,JJ+1) --> ZX
!             CAM0(II+1,JJ+2) --> XY
!             CAM0(II+2,JJ+2) --> YY
!             CAM0(II+3,JJ+2) --> ZY
!             CAM0(II+1,JJ+3) --> XZ
!             CAM0(II+2,JJ+3) --> YZ
!             CAM0(II+3,JJ+3) --> ZZ
        
        ENDDO   ! inner loop on dipoles
      ENDDO   ! outer loop on dipoles
!$OMP ENDDO

!$OMP END PARALLEL

      END SUBROUTINE SETUPMATRIX


      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

!********************************************************************
!--------------------------------------------------------------------*
!      Iversion of the matrix using a sequence of dyson equations    *
!                Returns the generalized Propagator K                *
!--------------------------------------------------------------------*
!*********************************************************************
      SUBROUTINE DYSONSEQ(ALAMBDA,SPACE, &
	          & XM,YM,ZM,  &
	          & CAP, &
	          & cn1,cn2,cn3, &
	          & CNORM, &
	          & NMAX, RETURNSUSC, &
	          & CAM0)
! 	          & NTHREADS, &
       
      !f2py threadsafe
      !$ use OMP_LIB
!       USE OMP_LIB
      USE PYGDMPRECISION, ONLY : dp
      implicit none

      
!**********************************************************************
!*      CONFIGURATION
!**********************************************************************
      ! Fixed allocations
      real(dp), parameter :: Pi=3.141592654_dp
      complex(dp), parameter :: C0=(0._dp,0._dp)
      complex(dp), parameter :: CUN=(1._dp,0._dp)
      complex(dp), parameter :: CIM=(0._dp,1._dp)
     
      ! parameters
      integer, intent(in) :: NMAX 
      integer, intent(in) :: RETURNSUSC
!       integer, intent(in) :: NTHREADS
      complex(dp),   intent(in) :: cn1,cn2,cn3
      real(dp), intent(in) :: ALAMBDA
      real(dp), intent(in) :: SPACE
      real(dp), intent(in) :: XM(NMAX),YM(NMAX),ZM(NMAX)
      complex(dp), intent(in) :: CAP(NMAX)  ! complex polarizabilities of each meshpoint
      complex(dp), intent(in) :: CNORM  ! mesh-type dependent normalization factor
      
      complex(dp), intent(out) :: CAM0(3*NMAX,3*NMAX)
     
      
      
!**********************************************************************
!*      DECLARATION
!**********************************************************************
!       integer :: NTHREADSUSE
      integer :: I,J,K,L,I1,J1,II,KK,JJ,LL,II1,JJ1
      real(dp) :: AK00
      
      real(dp) :: XI,YI,ZI,XJ,YJ,ZJ,ZS,XD,YD,ZD
      complex(dp) :: CHXX,CHXY,CHXZ,CHYX,CHYY,CHYZ,CHZX,CHZY,CHZZ
      complex(dp) :: SXX,SXY,SXZ,SYY,SYZ,SZZ
      complex(dp) :: CTXX,CTXY,CTXZ,CTYY,CTYZ,CTZZ
      complex(dp) :: CSXX,CSXY,CSXZ,CSYY,CSYZ,CSZZ
      real(dp) :: TXX,TXY,TXZ,TYY,TYZ,TZZ
      real(dp) :: TXX2,TXY2,TXZ2,TYY2,TYZ2,TZZ2
      
      complex(dp) :: CAXX(NMAX),CAYY(NMAX),CAZZ(NMAX),CAXY(NMAX)
      complex(dp) :: CAXZ(NMAX),CAYZ(NMAX),CAZX(NMAX),CAZY(NMAX),CAYX(NMAX)
      complex(dp) :: CVXX(NMAX),CVYY(NMAX),CVZZ(NMAX),CVXY(NMAX)
      complex(dp) :: CVXZ(NMAX),CVYX(NMAX),CVYZ(NMAX),CVZX(NMAX),CVZY(NMAX)
      complex(dp) :: CWWXX(NMAX),CWWYY(NMAX),CWWZZ(NMAX),CWWXY(NMAX)
      complex(dp) :: CWWXZ(NMAX),CWWYX(NMAX),CWWYZ(NMAX),CWWZX(NMAX),CWWZY(NMAX)
            
      complex(dp) :: CAAXX,CAAXY,CAAXZ,CAAYX,CAAYY,CAAYZ,CAAZX,CAAZY,CAAZZ
      complex(dp) :: CBXX,CBXY,CBXZ,CBYX,CBYY,CBYZ,CBZX,CBZY,CBZZ
      complex(dp) :: CXX,CXY,CXZ,CYX,CYY,CYZ,CZX,CZY,CZZ
      complex(dp) :: CDETERM,CDETER
      
     
     
      
      
      
!       IF (NTHREADS==-1) THEN
!         NTHREADSUSE = OMP_GET_MAX_THREADS()
!       ELSE
!         NTHREADSUSE = NTHREADS
!       ENDIF
!       CALL omp_set_num_threads(NTHREADSUSE)
 
 
!*********************************************************************
!*          INITIALIZATION
!*********************************************************************
 AK00=2._dp*Pi/ALAMBDA
     
  
!$OMP PARALLEL

!$OMP DO
!-----(1) Definition of the anisotropic polarizability of each -------
!---------------------discretized  cell ------------------------------
        DO J=1,NMAX
            CAXX(J)=CAP(J)
            CAYY(J)=CAP(J)
            CAZZ(J)=CAP(J)
            CAXY(J)=C0
            CAYX(J)=C0
            CAXZ(J)=C0
            CAZX(J)=C0
            CAZY(J)=C0
            CAYZ(J)=C0
        ENDDO
!$OMP ENDDO

!$OMP END PARALLEL


!---------Initialization of the matrix composed with ----------------
!-------------------site-to-site propagators-------------------------
!--------------------------AM0(I,J)----------------------------------
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(XJ,YJ,ZJ,JJ, XI,YI,ZI,II, XD,YD,ZD,ZS, &
!$OMP&              SXX,SYY,SZZ,SXY,SXZ,SYZ,  &
!$OMP&              CHXX,CHYY,CHZZ,CHXY,CHYX,CHYZ,CHZY,CHXZ,CHZX, &
!$OMP&              TXX,TYY,TZZ,TXY,TXZ,TYZ,TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2, &
!$OMP&              CTXX,CTYY,CTZZ,CTXY,CTXZ,CTYZ, &
!$OMP&              CSXX,CSYY,CSZZ,CSXY,CSXZ,CSYZ &
!$OMP&              )

!$OMP DO
      DO J=1,NMAX
        
        XJ=XM(J)
        YJ=YM(J)
        ZJ=ZM(J)
        JJ=3*(J-1)

        DO I=1,NMAX
            XI=XM(I)
            YI=YM(I)
            ZI=ZM(I)
            II=3*(I-1)

            IF(I.EQ.J)THEN
                ! --- "self"-term
                XD=0._dp
                YD=0._dp
                ZS=ZI+ZJ
                
                !! --- Surface propagator
                CALL PROPAS(XD,YD,ZS,SXX,SYY,SZZ,SXY,SXZ,SYZ,SPACE, &
                          & cn1,cn2,cn3) 
                
                !! Setup matrix G. Final matrix to inverse will be: (1 - G.alpha)
                CHXX = (SXX+CNORM)
                CHYY = (SYY+CNORM)
                CHZZ = (SZZ+CNORM)
                     
                CHXY =   SXY*CUN
                CHYX =   SXY*CUN
                CHYZ =   SYZ*CUN
                CHZY = (-1._dp*SYZ*CUN)
                CHXZ =   SXZ*CUN
                CHZX = (-1._dp*SXZ*CUN)
            ELSE
                ! --- Couple 2 different dipoles
                XD=XI-XJ
                YD=YI-YJ
                ZD=ZI-ZJ
                ZS=ZI+ZJ
                
                !! --- Vacuum propagator
                CALL PROPA0(AK00,XD,YD,ZD,TXX,TYY,TZZ,&
                          & TXY,TXZ,TYZ,TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2,&
                          & cn2)
                
                CTXX=TXX*CUN+TXX2*CIM
                CTYY=TYY*CUN+TYY2*CIM
                CTZZ=TZZ*CUN+TZZ2*CIM
                CTXY=TXY*CUN+TXY2*CIM
                CTXZ=TXZ*CUN+TXZ2*CIM
                CTYZ=TYZ*CUN+TYZ2*CIM
                
                !! --- Surface propagator
                CALL PROPAS(XD,YD,ZS,SXX,SYY,SZZ,SXY,SXZ,SYZ,SPACE,&
                          & cn1,cn2,cn3)
                
                CSXX=SXX*CUN
                CSYY=SYY*CUN
                CSZZ=SZZ*CUN
                CSXY=SXY*CUN
                CSXZ=SXZ*CUN
                CSYZ=SYZ*CUN
            !------------------------------------------------------------------
                CHXX = (CTXX+CSXX)
                CHYY = (CTYY+CSYY)
                CHZZ = (CTZZ+CSZZ)
                CHXY = (CTXY+CSXY)
                CHXZ = (CTXZ+CSXZ)
                CHYZ = (CTYZ+CSYZ)
                CHYX = (CTXY+CSXY)
                CHZY = (CTYZ-CSYZ)
                CHZX = (CTXZ-CSXZ)
                
            ENDIF
            
            CAM0(II+1,JJ+1)=CHXX
            CAM0(II+2,JJ+1)=CHYX
            CAM0(II+3,JJ+1)=CHZX
            
            CAM0(II+1,JJ+2)=CHXY
            CAM0(II+2,JJ+2)=CHYY
            CAM0(II+3,JJ+2)=CHZY
            
            CAM0(II+1,JJ+3)=CHXZ
            CAM0(II+2,JJ+3)=CHYZ
            CAM0(II+3,JJ+3)=CHZZ
            
        ENDDO   ! inner loop on dipoles
      ENDDO   ! outer loop on dipoles
!$OMP ENDDO

!$OMP END PARALLEL












!$OMP PARALLEL

!-------------- Matrix inversion via Dyson-sequence------------------*
!      We use in this code the concept of generalized propagator:   *
!cf. O.J.F. Martin, Ch. Girard, A. Dereux, Phys. Rev. Let. 74(1995)526
!           combined with a 3D sequence of Dyson's equations        *
!cf. Ch. Girard, O.J.F. Martin, A. Dereux,Phys. Rev. Let. 75(1995)3098
!********************************************************************
!
!
!----------SEQUENCE OF DYSON EQUATION---------------------------------
!
!--------------(1) Summation on the atoms-----------------------------
      DO K=1,NMAX
        KK=3*(K-1)
    !--------- Filling the 3x3 matrix A(K) -------------------------------
        CAAXX=1_dp-(CAM0(KK+1,KK+1)*CAXX(K)+CAM0(KK+1,KK+2)*CAYX(K)&
        &                            +CAM0(KK+1,KK+3)*CAZX(K))
        CAAYY=1_dp-(CAM0(KK+2,KK+1)*CAXY(K)+CAM0(KK+2,KK+2)*CAYY(K)&
        &                            +CAM0(KK+2,KK+3)*CAZY(K))
        CAAZZ=1_dp-(CAM0(KK+3,KK+1)*CAXZ(K)+CAM0(KK+3,KK+2)*CAYZ(K)&
        &                            +CAM0(KK+3,KK+3)*CAZZ(K))
        CAAXY=-(CAM0(KK+1,KK+1)*CAXY(K)+CAM0(KK+1,KK+2)*CAYY(K)&
        &                           +CAM0(KK+1,KK+3)*CAZY(K))
        CAAYX=-(CAM0(KK+2,KK+1)*CAXX(K)+CAM0(KK+2,KK+2)*CAYX(K)&
        &                           +CAM0(KK+2,KK+3)*CAZX(K))
        CAAXZ=-(CAM0(KK+1,KK+1)*CAXZ(K)+CAM0(KK+1,KK+2)*CAYZ(K)&
        &                           +CAM0(KK+1,KK+3)*CAZZ(K))
        CAAZX=-(CAM0(KK+3,KK+1)*CAXX(K)+CAM0(KK+3,KK+2)*CAYX(K)&
        &                           +CAM0(KK+3,KK+3)*CAZX(K))
        CAAYZ=-(CAM0(KK+2,KK+1)*CAXZ(K)+CAM0(KK+2,KK+2)*CAYZ(K)&
        &                           +CAM0(KK+2,KK+3)*CAZZ(K))
        CAAZY=-(CAM0(KK+3,KK+1)*CAXY(K)+CAM0(KK+3,KK+2)*CAYY(K)&
        &                           +CAM0(KK+3,KK+3)*CAZY(K))
     
!---------------------------------------------------------------------
!
!---------- Computing the inverse of AA(K) -> B(K) -------------------
!
!---------- Calculation of the determinant ---------------------------
        CDETERM =   CAAXX * (CAAYY*CAAZZ-CAAYZ*CAAZY) &
        &          - CAAYX * (CAAXY*CAAZZ-CAAZY*CAAXZ) &
        &          + CAAZX * (CAAXY*CAAYZ-CAAYY*CAAXZ)
    
        CDETER = CUN/CDETERM
      
      
!---------------------------------------------------------------------

!----------- Matrix inversion ----------------------------------------
        CBXX=CDETER*(CAAYY*CAAZZ-CAAZY*CAAYZ)
        CBXY=-CDETER*(CAAXY*CAAZZ-CAAZY*CAAXZ)
        CBXZ=CDETER*(CAAXY*CAAYZ-CAAYY*CAAXZ)
        CBYX=-CDETER*(CAAYX*CAAZZ-CAAZX*CAAYZ)
        CBYY=CDETER*(CAAXX*CAAZZ-CAAZX*CAAXZ)
        CBYZ=-CDETER*(CAAXX*CAAYZ-CAAYX*CAAXZ)
        CBZX=CDETER*(CAAYX*CAAZY-CAAZX*CAAYY)
        CBZY=-CDETER*(CAAXX*CAAZY-CAAZX*CAAXY)
        CBZZ=CDETER*(CAAXX*CAAYY-CAAYX*CAAXY) 
      
!---------------------------------------------------------------------
!
!------------ Matrix C(K) -> Alpha(K).B(K) ---------------------------
        CXX=CAXX(K)*CBXX+CAXY(K)*CBYX+CAXZ(K)*CBZX
        CYY=CAYX(K)*CBXY+CAYY(K)*CBYY+CAYZ(K)*CBZY
        CZZ=CAZX(K)*CBXZ+CAZY(K)*CBYZ+CAZZ(K)*CBZZ
        CXY=CAXX(K)*CBXY+CAXY(K)*CBYY+CAXZ(K)*CBZY
        CYX=CAYX(K)*CBXX+CAYY(K)*CBYX+CAYZ(K)*CBZX
        CXZ=CAXX(K)*CBXZ+CAXY(K)*CBYZ+CAXZ(K)*CBZZ
        CZX=CAZX(K)*CBXX+CAZY(K)*CBYX+CAZZ(K)*CBZX
        CYZ=CAYX(K)*CBXZ+CAYY(K)*CBYZ+CAYZ(K)*CBZZ
        CZY=CAZX(K)*CBXY+CAZY(K)*CBYY+CAZZ(K)*CBZY
      
!----------------------------------------------------------------------
!
      
!$OMP DO
      DO L=1,NMAX
        LL=3*(L-1)
    !----------------- Vector WW(L) ---------------------------------------
        CWWXX(L)=CAM0(LL+1,KK+1)*CXX+CAM0(LL+1,KK+2)*CYX&
        &       +CAM0(LL+1,KK+3)*CZX
        CWWYY(L)=CAM0(LL+2,KK+1)*CXY+CAM0(LL+2,KK+2)*CYY&
        &       +CAM0(LL+2,KK+3)*CZY
        CWWZZ(L)=CAM0(LL+3,KK+1)*CXZ+CAM0(LL+3,KK+2)*CYZ&
        &       +CAM0(LL+3,KK+3)*CZZ
        CWWXY(L)=CAM0(LL+1,KK+1)*CXY+CAM0(LL+1,KK+2)*CYY&
        &       +CAM0(LL+1,KK+3)*CZY
        CWWYX(L)=CAM0(LL+2,KK+1)*CXX+CAM0(LL+2,KK+2)*CYX&
        &       +CAM0(LL+2,KK+3)*CZX
        CWWXZ(L)=CAM0(LL+1,KK+1)*CXZ+CAM0(LL+1,KK+2)*CYZ&
        &       +CAM0(LL+1,KK+3)*CZZ
        CWWZX(L)=CAM0(LL+3,KK+1)*CXX+CAM0(LL+3,KK+2)*CYX& 
        &       +CAM0(LL+3,KK+3)*CZX
        CWWYZ(L)=CAM0(LL+2,KK+1)*CXZ+CAM0(LL+2,KK+2)*CYZ&
        &       +CAM0(LL+2,KK+3)*CZZ
        CWWZY(L)=CAM0(LL+3,KK+1)*CXY+CAM0(LL+3,KK+2)*CYY&
        &       +CAM0(LL+3,KK+3)*CZY
    !------------------ Vector V(L) ---------------------------------------- 
        CVXX(L)=CAM0(KK+1,LL+1)
        CVYY(L)=CAM0(KK+2,LL+2)
        CVZZ(L)=CAM0(KK+3,LL+3)
        CVXY(L)=CAM0(KK+1,LL+2)
        CVYX(L)=CAM0(KK+2,LL+1)
        CVXZ(L)=CAM0(KK+1,LL+3)
        CVZX(L)=CAM0(KK+3,LL+1)
        CVYZ(L)=CAM0(KK+2,LL+3)
        CVZY(L)=CAM0(KK+3,LL+2)

      ENDDO
!$OMP END DO
        
      
!-----------------------------------------------------------------------
!
!************************************************************************
!*       Solving the Dyson's equation inside the molecular              *
!*                           Structure                                  *
!************************************************************************
        DO J1=1,NMAX

!$OMP DO

            DO I1=1,NMAX
                II1=3*(I1-1)
                JJ1=3*(J1-1)
                !!!!!!!!!!!!! Order changed
                CAM0(II1+1,JJ1+1)=CAM0(II1+1,JJ1+1)+&
                &   CWWXX(I1)*CVXX(J1)+CWWXY(I1)*CVYX(J1)+CWWXZ(I1)*CVZX(J1)
                CAM0(II1+2,JJ1+1)=CAM0(II1+2,JJ1+1)+&
                &   CWWYX(I1)*CVXX(J1)+CWWYY(I1)*CVYX(J1)+CWWYZ(I1)*CVZX(J1)
                CAM0(II1+3,JJ1+1)=CAM0(II1+3,JJ1+1)+&
                &   CWWZX(I1)*CVXX(J1)+CWWZY(I1)*CVYX(J1)+CWWZZ(I1)*CVZX(J1)
                CAM0(II1+1,JJ1+2)=CAM0(II1+1,JJ1+2)+&
                &   CWWXX(I1)*CVXY(J1)+CWWXY(I1)*CVYY(J1)+CWWXZ(I1)*CVZY(J1)
                CAM0(II1+2,JJ1+2)=CAM0(II1+2,JJ1+2)+&
                &   CWWYX(I1)*CVXY(J1)+CWWYY(I1)*CVYY(J1)+CWWYZ(I1)*CVZY(J1)
                CAM0(II1+3,JJ1+2)=CAM0(II1+3,JJ1+2)+&
                &   CWWZX(I1)*CVXY(J1)+CWWZY(I1)*CVYY(J1)+CWWZZ(I1)*CVZY(J1)
                CAM0(II1+1,JJ1+3)=CAM0(II1+1,JJ1+3)+&
                &   CWWXX(I1)*CVXZ(J1)+CWWXY(I1)*CVYZ(J1)+CWWXZ(I1)*CVZZ(J1)
                CAM0(II1+2,JJ1+3)=CAM0(II1+2,JJ1+3)+&
                &   CWWYX(I1)*CVXZ(J1)+CWWYY(I1)*CVYZ(J1)+CWWYZ(I1)*CVZZ(J1)
                CAM0(II1+3,JJ1+3)=CAM0(II1+3,JJ1+3)+&
                &   CWWZX(I1)*CVXZ(J1)+CWWZY(I1)*CVYZ(J1)+CWWZZ(I1)*CVZZ(J1)

            ENDDO
!$OMP END DO

        ENDDO
      ENDDO
      
!$OMP END PARALLEL
!************************************************************************
!*          End of the Dyson's sequence  --> now CAM0 is field suscept. *
!************************************************************************
!*                                                                      *
!*                                                                      *
!******************** Calc. Generalized Propagator **********************
    IF (RETURNSUSC /= 1) THEN
!$OMP PARALLEL      
      DO K=1,NMAX
!$OMP DO      
        DO L=1,NMAX
            LL=3*(L-1)
            KK=3*(K-1) 
            IF(L.EQ.K)THEN
                CAM0(LL+1,KK+1)=1+CAM0(LL+1,KK+1)*CAXX(K)+CAM0(LL+1,KK+2)*CAYX(K)&
                &       +CAM0(LL+1,KK+3)*CAZX(K)
                CAM0(LL+2,KK+1)=CAM0(LL+2,KK+1)*CAXX(K)+CAM0(LL+2,KK+2)*CAYX(K)&
                &       +CAM0(LL+2,KK+3)*CAZX(K)
                CAM0(LL+3,KK+1)=CAM0(LL+3,KK+1)*CAXX(K)+CAM0(LL+3,KK+2)*CAYX(K)&
                &       +CAM0(LL+3,KK+3)*CAZX(K)
                CAM0(LL+1,KK+2)=CAM0(LL+1,KK+1)*CAXY(K)+CAM0(LL+1,KK+2)*CAYY(K)&
                &       +CAM0(LL+1,KK+3)*CAZY(K)
                
                CAM0(LL+2,KK+2)=1+CAM0(LL+2,KK+1)*CAXY(K)+CAM0(LL+2,KK+2)*CAYY(K)&
                &       +CAM0(LL+2,KK+3)*CAZY(K)
                CAM0(LL+3,KK+2)=CAM0(LL+3,KK+1)*CAXY(K)+CAM0(LL+3,KK+2)*CAYY(K)&
                &       +CAM0(LL+3,KK+3)*CAZY(K)
                CAM0(LL+1,KK+3)=CAM0(LL+1,KK+1)*CAXZ(K)+CAM0(LL+1,KK+2)*CAYZ(K)&
                &       +CAM0(LL+1,KK+3)*CAZZ(K)
                CAM0(LL+2,KK+3)=CAM0(LL+2,KK+1)*CAXZ(K)+CAM0(LL+2,KK+2)*CAYZ(K)&
                &       +CAM0(LL+2,KK+3)*CAZZ(K)
                
                CAM0(LL+3,KK+3)=1+CAM0(LL+3,KK+1)*CAXZ(K)+CAM0(LL+3,KK+2)*CAYZ(K)&
                &       +CAM0(LL+3,KK+3)*CAZZ(K)
            ELSE
                CAM0(LL+1,KK+1)=CAM0(LL+1,KK+1)*CAXX(K)+CAM0(LL+1,KK+2)*CAYX(K)&
                &       +CAM0(LL+1,KK+3)*CAZX(K)
                CAM0(LL+2,KK+1)=CAM0(LL+2,KK+1)*CAXX(K)+CAM0(LL+2,KK+2)*CAYX(K)&
                &       +CAM0(LL+2,KK+3)*CAZX(K)
                CAM0(LL+3,KK+1)=CAM0(LL+3,KK+1)*CAXX(K)+CAM0(LL+3,KK+2)*CAYX(K)&
                &       +CAM0(LL+3,KK+3)*CAZX(K)
                
                CAM0(LL+1,KK+2)=CAM0(LL+1,KK+1)*CAXY(K)+CAM0(LL+1,KK+2)*CAYY(K)&
                &       +CAM0(LL+1,KK+3)*CAZY(K)
                CAM0(LL+2,KK+2)=CAM0(LL+2,KK+1)*CAXY(K)+CAM0(LL+2,KK+2)*CAYY(K)&
                &       +CAM0(LL+2,KK+3)*CAZY(K)
                CAM0(LL+3,KK+2)=CAM0(LL+3,KK+1)*CAXY(K)+CAM0(LL+3,KK+2)*CAYY(K)&
                &       +CAM0(LL+3,KK+3)*CAZY(K)
                
                CAM0(LL+1,KK+3)=CAM0(LL+1,KK+1)*CAXZ(K)+CAM0(LL+1,KK+2)*CAYZ(K)&
                &       +CAM0(LL+1,KK+3)*CAZZ(K)
                CAM0(LL+2,KK+3)=CAM0(LL+2,KK+1)*CAXZ(K)+CAM0(LL+2,KK+2)*CAYZ(K)&
                &       +CAM0(LL+2,KK+3)*CAZZ(K)
                CAM0(LL+3,KK+3)=CAM0(LL+3,KK+1)*CAXZ(K)+CAM0(LL+3,KK+2)*CAYZ(K)&
                &       +CAM0(LL+3,KK+3)*CAZZ(K)
            ENDIF
        ENDDO
!$OMP END DO      
      ENDDO
!$OMP END PARALLEL    
    ENDIF

!**************** REMARK: At this stage the matrix CAM0 is ***************
!***************** filled with the generalized propagator ***************

      END SUBROUTINE DYSONSEQ



      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
