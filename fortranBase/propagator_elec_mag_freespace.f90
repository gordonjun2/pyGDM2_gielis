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
! Author: Christian Girard 
! CEMES/CNRS 
!
! Date: January 1997
!
! ***********
! Modified 2016 by P.R. Wiecha
! CEMES/CNRS
!   - variable precision
!
!*********************************************************************


!C*****Christian GIRARD, 24 Janvier 1997***************************
      SUBROUTINE ELECMAG0(AK0,XC,YC,ZC,EBXY1,EBXZ1,EBYZ1, &
                                               & EBXY2,EBXZ2,EBYZ2)
     
      USE PYGDMPRECISION, ONLY : dp
      implicit none
      
      complex(dp), parameter :: c0=(0._dp,0._dp)
      complex(dp), parameter :: cun=(1._dp,0._dp)
      complex(dp), parameter :: cim=(0._dp,1._dp)
      real(dp), parameter :: Pi=3.141592654_dp
      
      
      real(dp) :: AK0,AK02
      real(dp) :: XC,YC,ZC
      real(dp) :: EBXY1,EBXZ1,EBYZ1
      real(dp) :: EBXY2,EBXZ2,EBYZ2
      complex(dp) :: CK0
      real(dp) :: T2XY,T2YZ,T2XZ
      real(dp) :: T3XY,T3YZ,T3XZ
      real(dp) :: RMO
      complex(dp) :: CFEXP
      complex(dp) :: CTXY,CTYZ,CTXZ
      
      
      AK02=AK0**2
      CK0=CIM*AK0
!C*****************************************************************
!C******CONSTRUCTION DES PROPAGATEURS RETARDES  *******************
!C*****************************************************************
      RMO=(XC*XC+YC*YC+ZC*ZC)
!C-----------------------------------------------------------------
      T2XY=-ZC/RMO
      T3XY=-ZC/RMO**1.5
!C-----------------------------------------------------------------
      T2XZ=YC/RMO
      T3XZ=YC/RMO**1.5
!C-----------------------------------------------------------------
      T2YZ=-XC/RMO
      T3YZ=-XC/RMO**1.5
!C-----------------------------------------------------------------
      CFEXP=EXP(CIM*AK0*SQRT(RMO))
      
      CTXY=CFEXP*(CK0*T3XY+AK02*T2XY)
      CTXZ=CFEXP*(CK0*T3XZ+AK02*T2XZ)
      CTYZ=CFEXP*(CK0*T3YZ+AK02*T2YZ)
      
      EBXY1=REAL(CTXY)
      EBXZ1=REAL(CTXZ)
      EBYZ1=REAL(CTYZ)
      EBXY2=AIMAG(CTXY)
      EBXZ2=AIMAG(CTXZ)
      EBYZ2=AIMAG(CTYZ)
      END





!! --- helper-routine to get full S^EB tensor (complex)
!! --- input parameters follow the convention: S^EB (R_OBS, R_C)
SUBROUTINE ELECMAG0_FULL(AK0,XC,YC,ZC, XOBS,YOBS,ZOBS, &
              & CEBXX, CEBXY, CEBXZ, &
              & CEBYX, CEBYY, CEBYZ, &
              & CEBZX, CEBZY, CEBZZ, &
              & cn1,cn2,cn3)
     
      USE PYGDMPRECISION, ONLY : dp
      implicit none
      
      complex(dp), parameter :: c0=(0._dp,0._dp)
      complex(dp), parameter :: cun=(1._dp,0._dp)
      complex(dp), parameter :: cim=(0._dp,1._dp)
      real(dp), parameter :: Pi=3.141592654_dp
      
      
      real(dp), intent(in) :: AK0
      real(dp), intent(in) :: XC,YC,ZC
      real(dp), intent(in) :: XOBS,YOBS,ZOBS
      real(dp) :: XDDD,YDDD,ZDDD
      
      real(dp) :: EBXY1,EBXZ1,EBYZ1
      real(dp) :: EBXY2,EBXZ2,EBYZ2
!       real(dp) :: EBXX1,EBXY1,EBXZ1
!       real(dp) :: EBYX1,EBYY1,EBYZ1
!       real(dp) :: EBZX1,EBZY1,EBZZ1
!       real(dp) :: EBXX2,EBXY2,EBXZ2
!       real(dp) :: EBYX2,EBYY2,EBYZ2
!       real(dp) :: EBZX2,EBZY2,EBZZ2
      
      complex(dp), intent(in) :: cn1,cn2,cn3
      
!       complex(dp) :: CTXY,CTYZ,CTXZ
      
      complex(dp), intent(out) :: CEBXX, CEBYY, CEBZZ
      complex(dp), intent(out) :: CEBXY, CEBYX, CEBXZ
      complex(dp), intent(out) :: CEBZX, CEBYZ, CEBZY
      
!       if ( (cn1.ne.cn2).or.(cn1.ne.cn3) ) then
!        write (*,*) "Error. Only defined for homogeneous media! please use according refractive-indices."
!       end if
      
      XDDD = XOBS - XC
      YDDD = YOBS - YC
      ZDDD = ZOBS - ZC
      
        CALL ELECMAG0(AK0,XDDD,YDDD,ZDDD,&
            & EBXY1,EBXZ1,EBYZ1,&
            & EBXY2,EBXZ2,EBYZ2)
        CEBXX=C0
        CEBYY=C0
        CEBZZ=C0
        CEBXY=EBXY1*CUN+EBXY2*CIM
        CEBYX=-CEBXY
        CEBXZ=EBXZ1*CUN+EBXZ2*CIM
        CEBZX=-CEBXZ
        CEBYZ=EBYZ1*CUN+EBYZ2*CIM
        CEBZY=-CEBYZ
      END 
        
        
        

!! --- helper-routine to get full S^EB tensor (complex)
!! --- input parameters follow the convention: S^BB (R_OBS, R_C)
SUBROUTINE MAGMAG0_FULL(AK0, XOBS,YOBS,ZOBS, XC,YC,ZC, &
              & CBBXX, CBBXY, CBBXZ, &
              & CBBYX, CBBYY, CBBYZ, &
              & CBBZX, CBBZY, CBBZZ, &
              & cn1,cn2,cn3)
     
      USE PYGDMPRECISION, ONLY : dp
      implicit none
      
      complex(dp), parameter :: c0=(0._dp,0._dp)
      complex(dp), parameter :: cun=(1._dp,0._dp)
      complex(dp), parameter :: cim=(0._dp,1._dp)
      real(dp), parameter :: Pi=3.141592654_dp
      
      real(dp), intent(in) :: AK0
      real(dp), intent(in) :: XC,YC,ZC
       real(dp), intent(in) :: XOBS,YOBS,ZOBS
      real(dp) :: ALAMBDA, SPACING
      
      complex(dp), intent(in) :: cn1,cn2,cn3
      
      complex(dp), intent(out) :: CBBXX, CBBYY, CBBZZ
      complex(dp), intent(out) :: CBBXY, CBBYX, CBBXZ
      complex(dp), intent(out) :: CBBZX, CBBYZ, CBBZY
      
      ALAMBDA = 2._dp*Pi/AK0
      SPACING = 5000
      
        CALL PROPAGATOR123(ALAMBDA,&
                & XOBS,YOBS,ZOBS,XC,YC,ZC, SPACING,&
                & CBBXX,CBBXY,CBBXZ,CBBYX,CBBYY,CBBYZ,CBBZX,CBBZY,CBBZZ, &
                & cn1,cn2,cn3)
        
      END 
!!    -----------------------------------------------------------------------



