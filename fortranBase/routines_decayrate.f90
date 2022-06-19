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
!!*********************************************************************
! Author: Peter Wiecha
! CEMES/CNRS 
!
! Date: July 2016
!
!*********************************************************************







!! --- routine to get full S_P^EE / S_P^BB tensor (complex)
!! --- requires pre-calculated generalized propagator CAM0 (usually called "K")
      SUBROUTINE PROPA_SP_DECAY( &
                    & STEP,NMAX,XM,YM,ZM, N1,N2,N3, &
                    & ALAMBDA, &
                    & NPTSIN, XMAP,YMAP,ZMAP, &
                    & MAGNETIC, &
                    & CAM0, ALPHA, SP)
!                     & NTHREADS, &
      
      !f2py threadsafe
      !$ use OMP_LIB
!       use OMP_LIB
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
      
      !***************************************
      !*      Structure properties  
      !***************************************
      !!! ATTENTION: Structure must be shifted to    Z(min) >= 0.5*step !!!
      real(dp), intent(in) :: STEP
      real(dp) :: VCELL
      
      integer, intent(in) :: NMAX
      real(dp), intent(in) :: XM(NMAX), YM(NMAX), ZM(NMAX)
      real(dp), intent(in) :: N1,N2,N3    ! RefIndex: SUBSTRATE(1), STRUCTURE SOURROUNDING(2), TOP MEDIUM(3)
      complex(dp), intent(in) :: ALPHA
      complex(dp) :: CN1,CN2,CN3
      
      
      !***************************************
      !*      Simulation
      !***************************************
!       integer, intent(in) :: NTHREADS
      
      ! wavelength (nm)
      real(dp), intent(in)  :: ALAMBDA
      real(dp)  :: AK0
      
      ! simulation type, 1=gamma_m, else=gamma_e
      integer, intent(in)  :: MAGNETIC
      
      
      !***************************************
      !*      NF-"Carte" Definition
      !***************************************
      integer, intent(in) :: NPTSIN
      real(dp), intent(in)    :: XMAP(NPTSIN),YMAP(NPTSIN),ZMAP(NPTSIN)
      
      
      
      !***************************************
      !*      Propagators
      !***************************************
      complex(dp), intent(in) :: CAM0(3*NMAX, 3*NMAX)
      complex(dp), intent(out) :: SP(NPTSIN, 9)

      
      
      
      
      !***************************************
      !*      OTHER DECLARATIONS
      !***************************************
!       integer :: NTHREADSUSE
      
      !* ----------- Temporary Propagators
      complex(dp) :: Q0XX,Q0YY,Q0ZZ, Q0XY,Q0YX,Q0XZ, Q0ZX,Q0YZ,Q0ZY
      complex(dp) :: S0XX,S0YY,S0ZZ, S0XY,S0YX,S0XZ, S0ZX,S0YZ,S0ZY
      complex(dp) :: KijXX,KijYY,KijZZ, KijXY,KijYX,KijXZ, KijZX,KijYZ,KijZY
      complex(dp) :: StmpXX,StmpYY,StmpZZ, StmpXY,StmpYX,StmpXZ, StmpZX,StmpYZ,StmpZY
      complex(dp) :: SPXX,SPYY,SPZZ, SPXY,SPYX,SPXZ, SPZX,SPYZ,SPZY
      
      real(dp) :: XObs,YObs,ZObs,  Xi,Yi,Zi,  Xj,Yj,Zj
      real(dp) :: XD,YD,ZD, XDDD,YDDD,ZDDD, SBB, SHORTER
      
      integer :: I,J, IMAP, INSIDE
      
      
      
      
      
!**********************************************************************
!*      Configure
!**********************************************************************
!       IF (NTHREADS==-1) THEN
!         NTHREADSUSE = OMP_GET_MAX_THREADS()
!       ELSE
!         NTHREADSUSE = NTHREADS
!       ENDIF
!       CALL omp_set_num_threads(NTHREADSUSE)
      
      VCELL = STEP**3
!* Sourrounding refindex
      CN1 = CUN*N1
      CN2 = CUN*N2
      CN3 = CUN*N3
      
!* Wavelength/Polarization
      AK0=2._dp*Pi/ALAMBDA
      
      
      !* --- Loop Cartography-points 
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(ALAMBDA,AK0, VCELL, CN1,CN2,CN3), &
!$OMP& SHARED(NMAX, NPTSIN, MAGNETIC), &
!$OMP& SHARED(XM,YM,ZM, XMAP,YMAP,ZMAP), &
!$OMP& SHARED(SP, CAM0)
!$OMP DO
      DO IMAP=1,NPTSIN
        XObs = XMAP(IMAP)
        YObs = YMAP(IMAP)
        ZObs = ZMAP(IMAP)
        
        !* --- initialize final matrix SP
        SP(IMAP, 1) = C0
        SP(IMAP, 2) = C0
        SP(IMAP, 3) = C0
        SP(IMAP, 4) = C0
        SP(IMAP, 5) = C0
        SP(IMAP, 6) = C0
        SP(IMAP, 7) = C0
        SP(IMAP, 8) = C0
        SP(IMAP, 9) = C0
        
        !* --- Check if inside or outside of structure
        SBB=STEP*0.999
      INSIDE = 0
      DO J=1,NMAX
        XD = XM(J)
        YD = YM(J)
        ZD = ZM(J)
        XDDD = XObs-XD
        YDDD = YObs-YD
        ZDDD = ZObs-ZD
        SHORTER=SQRT(XDDD**2+YDDD**2+ZDDD**2) 
        IF (SHORTER.LE.SBB) THEN
            !* --- inside structure: get closest meshpoint
            SBB=SHORTER
            INSIDE = J
        ENDIF
      ENDDO        
        
!!!*** =======================================
!!!*** INSIDE STRUCTURE
!!!*** =======================================
        IF (INSIDE /= 0) THEN
            !* --- take simply the diagonal tensor on generalized propagator for cell "inside"
            I = INSIDE
            J = INSIDE
            IF (MAGNETIC==1) THEN
                SP(IMAP, 1) = 0
                SP(IMAP, 4) = 0
                SP(IMAP, 7) = 0
                
                SP(IMAP, 2) = 0
                SP(IMAP, 5) = 0
                SP(IMAP, 8) = 0
                
                SP(IMAP, 3) = 0
                SP(IMAP, 6) = 0
                SP(IMAP, 9) = 0
            ELSE
                SP(IMAP, 1) = (CAM0(3*(I-1)+1, 3*(J-1)+1) - 1) / alpha
                SP(IMAP, 4) = CAM0(3*(I-1)+2, 3*(J-1)+1) / alpha
                SP(IMAP, 7) = CAM0(3*(I-1)+3, 3*(J-1)+1) / alpha
                
                SP(IMAP, 2) = CAM0(3*(I-1)+1, 3*(J-1)+2) / alpha
                SP(IMAP, 5) = (CAM0(3*(I-1)+2, 3*(J-1)+2) - 1) / alpha
                SP(IMAP, 8) = CAM0(3*(I-1)+3, 3*(J-1)+2) / alpha
                
                SP(IMAP, 3) = CAM0(3*(I-1)+1, 3*(J-1)+3) / alpha
                SP(IMAP, 6) = CAM0(3*(I-1)+2, 3*(J-1)+3) / alpha
                SP(IMAP, 9) = (CAM0(3*(I-1)+3, 3*(J-1)+3) - 1) / alpha
            ENDIF
            
!!!*** =======================================
!!!*** OUTSIDE STRUCTURE
!!!*** =======================================
        ELSEIF (INSIDE == 0) THEN
        !* --- Double-integral over structure
        DO I=1,NMAX
          Xi = XM(I)
          Yi = YM(I)
          Zi = ZM(I)
          DO J=1,NMAX
            Xj = XM(J)
            Yj = YM(J)
            Zj = ZM(J)
!! ELECMAG0_FULL: B-dipole; MAGMAG0_FULL: E-dipole
            IF (MAGNETIC==1) THEN
              CALL ELECMAG0_FULL(AK0, XObs,YObs,ZObs, Xi,Yi,Zi, &
                            & Q0XX, Q0XY, Q0XZ, &
                            & Q0YX, Q0YY, Q0YZ, &
                            & Q0ZX, Q0ZY, Q0ZZ, &
                            & cn1,cn2,cn3)
              Q0XX = -1*Q0XX
              Q0XY = -1*Q0XY
              Q0XZ = -1*Q0XZ
              Q0YX = -1*Q0YX
              Q0YY = -1*Q0YY
              Q0YZ = -1*Q0YZ
              Q0ZX = -1*Q0ZX
              Q0ZY = -1*Q0ZY
              Q0ZZ = -1*Q0ZZ
              CALL ELECMAG0_FULL(AK0, Xj,Yj,Zj, XObs,YObs,ZObs, &
                            & S0XX, S0XY, S0XZ, &
                            & S0YX, S0YY, S0YZ, &
                            & S0ZX, S0ZY, S0ZZ, &
                            & cn1,cn2,cn3)
            ELSE
              CALL MAGMAG0_FULL(AK0, XObs,YObs,ZObs, Xi,Yi,Zi, &
                            & Q0XX, Q0XY, Q0XZ, &
                            & Q0YX, Q0YY, Q0YZ, &
                            & Q0ZX, Q0ZY, Q0ZZ, &
                            & cn1,cn2,cn3)
              CALL MAGMAG0_FULL(AK0, Xj,Yj,Zj, XObs,YObs,ZObs, &
                            & S0XX, S0XY, S0XZ, &
                            & S0YX, S0YY, S0YZ, &
                            & S0ZX, S0ZY, S0ZZ, &
                            & cn1,cn2,cn3)
            ENDIF
            
            KijXX = CAM0(3*(I-1)+1, 3*(J-1)+1)
            KijYX = CAM0(3*(I-1)+2, 3*(J-1)+1)
            KijZX = CAM0(3*(I-1)+3, 3*(J-1)+1)
            KijXY = CAM0(3*(I-1)+1, 3*(J-1)+2)
            KijYY = CAM0(3*(I-1)+2, 3*(J-1)+2)
            KijZY = CAM0(3*(I-1)+3, 3*(J-1)+2)
            KijXZ = CAM0(3*(I-1)+1, 3*(J-1)+3)
            KijYZ = CAM0(3*(I-1)+2, 3*(J-1)+3)
            KijZZ = CAM0(3*(I-1)+3, 3*(J-1)+3)
            
            
            !* --- Calculate Stmp = K.S0
            StmpXX = KijXX*S0XX + KijXY*S0YX + KijXZ*S0ZX
            StmpXY = KijXX*S0XY + KijXY*S0YY + KijXZ*S0ZY
            StmpXZ = KijXX*S0XZ + KijXY*S0YZ + KijXZ*S0ZZ
            
            StmpYX = KijYX*S0XX + KijYY*S0YX + KijYZ*S0ZX
            StmpYY = KijYX*S0XY + KijYY*S0YY + KijYZ*S0ZY
            StmpYZ = KijYX*S0XZ + KijYY*S0YZ + KijYZ*S0ZZ
            
            StmpZX = KijZX*S0XX + KijZY*S0YX + KijZZ*S0ZX
            StmpZY = KijZX*S0XY + KijZY*S0YY + KijZZ*S0ZY
            StmpZZ = KijZX*S0XZ + KijZY*S0YZ + KijZZ*S0ZZ
            
            
            !* --- Calculate SP = Q0.Stmp
            SPXX = Q0XX*StmpXX + Q0XY*StmpYX + Q0XZ*StmpZX
            SPXY = Q0XX*StmpXY + Q0XY*StmpYY + Q0XZ*StmpZY
            SPXZ = Q0XX*StmpXZ + Q0XY*StmpYZ + Q0XZ*StmpZZ
            
            SPYX = Q0YX*StmpXX + Q0YY*StmpYX + Q0YZ*StmpZX
            SPYY = Q0YX*StmpXY + Q0YY*StmpYY + Q0YZ*StmpZY
            SPYZ = Q0YX*StmpXZ + Q0YY*StmpYZ + Q0YZ*StmpZZ
            
            SPZX = Q0ZX*StmpXX + Q0ZY*StmpYX + Q0ZZ*StmpZX
            SPZY = Q0ZX*StmpXY + Q0ZY*StmpYY + Q0ZZ*StmpZY
            SPZZ = Q0ZX*StmpXZ + Q0ZY*StmpYZ + Q0ZZ*StmpZZ
            
            
            !* --- construct final matrix
            SP(IMAP, 1) = SP(IMAP, 1) + SPXX
            SP(IMAP, 2) = SP(IMAP, 2) + SPXY
            SP(IMAP, 3) = SP(IMAP, 3) + SPXZ
            SP(IMAP, 4) = SP(IMAP, 4) + SPYX
            SP(IMAP, 5) = SP(IMAP, 5) + SPYY
            SP(IMAP, 6) = SP(IMAP, 6) + SPYZ
            SP(IMAP, 7) = SP(IMAP, 7) + SPZX
            SP(IMAP, 8) = SP(IMAP, 8) + SPZY
            SP(IMAP, 9) = SP(IMAP, 9) + SPZZ
            
          ENDDO
        ENDDO
!         write (*,*) "outside:", SP(IMAP,1)
        ENDIF
        
        !* --- normalization
        SP(IMAP, 1) = SP(IMAP, 1) * VCELL
        SP(IMAP, 2) = SP(IMAP, 2) * VCELL
        SP(IMAP, 3) = SP(IMAP, 3) * VCELL
        SP(IMAP, 4) = SP(IMAP, 4) * VCELL
        SP(IMAP, 5) = SP(IMAP, 5) * VCELL
        SP(IMAP, 6) = SP(IMAP, 6) * VCELL
        SP(IMAP, 7) = SP(IMAP, 7) * VCELL
        SP(IMAP, 8) = SP(IMAP, 8) * VCELL
        SP(IMAP, 9) = SP(IMAP, 9) * VCELL
      ENDDO
!$OMP END DO
!$OMP END PARALLEL 

      END



