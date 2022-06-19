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
!C*****Christian GIRARD, 14 Fevrier 1995****************************
!C*****Modified 9 April 1997 ***************************************
!
!   Free-space propagator (T1, T2, T3)
!
!C*************************************************************
     SUBROUTINE PROPA0(AK0,XC,YC,ZC,TXX1,TYY1,TZZ1, &
         & TXY1,TXZ1,TYZ1,TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2, &
         & cn2)
     
     USE PYGDMPRECISION, ONLY : dp
     IMPLICIT real(dp) (A-Z)
     
     
      complex(dp), parameter :: C0=(0._dp,0._dp)
      complex(dp), parameter :: CUN=(1._dp,0._dp)
      complex(dp), parameter :: CIM=(0._dp,1._dp)
      real(dp), parameter :: PI=3.141592654_dp
     
      real(dp), intent(in) :: AK0
      real(dp), intent(in) :: XC,YC,ZC
      real(dp), intent(out) :: TXX1,TYY1,TZZ1   ! real part
      real(dp), intent(out) :: TXY1,TXZ1,TYZ1   ! real part
      real(dp), intent(out) :: TXX2,TYY2,TZZ2   ! imag part
      real(dp), intent(out) :: TXY2,TXZ2,TYZ2   ! imag part
      real(dp) :: T1XX,T1XY,T1XZ,T1YY,T1YZ,T1ZZ
      real(dp) :: T2XX,T2XY,T2XZ,T2YY,T2YZ,T2ZZ
      real(dp) :: T3XX,T3XY,T3XZ,T3YY,T3YZ,T3ZZ   
      complex(dp), intent(in) :: cn2
      complex(dp) :: ceps2,CAK22,CK0
      real(dp) :: RMO
      complex(dp) :: CFEXP,CTXX,CTXY,CTXZ,CTYY,CTYZ,CTZZ
      
         
      
!C*****************************************************************
!C******CONSTRUCTION DES PROPAGATEURS RETARDES  ****************
!C*****************************************************************
      T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
 
     
      ceps2=cn2*cn2
      CAK22=AK0*AK0*ceps2
      CK0=-CIM*AK0*cn2
      RMO=SQRT(XC**2+YC**2+ZC**2)
      CFEXP=EXP(CIM*AK0*cn2*RMO)
      CTXX=CFEXP*(T1XX+CK0*T2XX-CAK22*T3XX)/ceps2
      CTYY=CFEXP*(T1YY+CK0*T2YY-CAK22*T3YY)/ceps2
      CTZZ=CFEXP*(T1ZZ+CK0*T2ZZ-CAK22*T3ZZ)/ceps2
      CTXY=CFEXP*(T1XY+CK0*T2XY-CAK22*T3XY)/ceps2
      CTXZ=CFEXP*(T1XZ+CK0*T2XZ-CAK22*T3XZ)/ceps2
      CTYZ=CFEXP*(T1YZ+CK0*T2YZ-CAK22*T3YZ)/ceps2
      TXX1=REAL(CTXX)
      TYY1=REAL(CTYY)
      TZZ1=REAL(CTZZ)
      TXY1=REAL(CTXY)
      TXZ1=REAL(CTXZ)
      TYZ1=REAL(CTYZ)
      TXX2=AIMAG(CTXX)
      TYY2=AIMAG(CTYY)
      TZZ2=AIMAG(CTZZ)
      TXY2=AIMAG(CTXY)
      TXZ2=AIMAG(CTXZ)
      TYZ2=AIMAG(CTYZ)
      END













!C*************************************************************
!C*****Ch Girard, 15 February 1996 ****************************
!C*************************************************************
!
!   Non-retarded surface propagator
!
!C*************************************************************
      SUBROUTINE PROPAS(X,Y,Z,SXX,SYY,SZZ,SXY,SXZ,SYZ,D,&
     &cn1,cn2,cn3)
     
     USE PYGDMPRECISION, ONLY : dp
     IMPLICIT real(dp) (A-Z)
     
     
      complex(dp), parameter :: C0=(0._dp,0._dp)
      complex(dp), parameter :: CUN=(1._dp,0._dp)
      complex(dp), parameter :: CIM=(0._dp,1._dp)
      complex(dp), intent(in) :: cn1,cn2,cn3      
      real(dp), intent(in) :: X,Y,Z,D
     complex(dp) :: ceps1,ceps2,ceps3
     complex(dp) :: cdelta12,cdelta23
     complex(dp) :: SXX12,SYY12,SZZ12,SXY12,SXZ12,SYZ12
     complex(dp) :: SXX23,SYY23,SZZ23,SXY23,SXZ23,SYZ23
     complex(dp), intent(out) :: SXX,SYY,SZZ,SXY,SXZ,SYZ
!************* Optical indexes ********************************
!*
!************* Real case **************************************
     
      ceps1=cn1**2
      ceps2=cn2**2
      ceps3=cn3**2
      cdelta12=(ceps1-ceps2)/(ceps1+ceps2)
      cdelta23=(ceps3-ceps2)/(ceps3+ceps2)
!       write (*,*) ceps1, cdelta12,cdelta23
!*
!************* Interface: (1,2) *******************************
      SXX12=cdelta12*(Z*Z+Y*Y-2_dp*X*X)/(X*X+Y*Y+Z*Z)**2.5
!C*************************************
      SYY12=cdelta12*(Z*Z+X*X-2_dp*Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!C*************************************
      SZZ12=cdelta12*(2_dp*Z*Z-X*X-Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!C**************************************
      SXY12=cdelta12*(-3_dp*X*Y)/(X*X+Y*Y+Z*Z)**2.5
!C***************************************
      SXZ12=cdelta12*(3_dp*X*Z)/(X*X+Y*Y+Z*Z)**2.5
!C***************************************
      SYZ12=cdelta12*(3_dp*Y*Z)/(X*X+Y*Y+Z*Z)**2.5
!**************************************************************
!*
!************* Interface: (2,3) *******************************
      GZ=Z-2._dp*D
      SXX23=cdelta23*(GZ*GZ+Y*Y-2_dp*X*X)/(X*X+Y*Y+GZ*GZ)**2.5
!C*************************************
      SYY23=cdelta23*(GZ*GZ+X*X-2_dp*Y*Y)/(X*X+Y*Y+GZ*GZ)**2.5
!C*************************************
      SZZ23=cdelta23*(2_dp*GZ*GZ-X*X-Y*Y)/(X*X+Y*Y+GZ*GZ)**2.5
!C**************************************
      SXY23=cdelta23*(-3_dp*X*Y)/(X*X+Y*Y+GZ*GZ)**2.5
!C***************************************
      SXZ23=cdelta23*(3_dp*X*GZ)/(X*X+Y*Y+GZ*GZ)**2.5
!C***************************************
      SYZ23=cdelta23*(3_dp*Y*GZ)/(X*X+Y*Y+GZ*GZ)**2.5
!*
!**************************************************************
      SXX=SXX12+SXX23
      SYY=SYY12+SYY23
      SZZ=SZZ12+SZZ23
      SXY=SXY12+SXY23
      SXZ=SXZ12+SXZ23
      SYZ=SYZ12+SYZ23
      END 
 
  
 













!C******************************************************************
!C****************** propagator_123.f ******************************
!C******************************************************************
!C****************** Author: Ch. Girard ****************************
!C****************** Date: 9 April 1997 ****************************
!C****************** Based on the retarded image theory ************
!C******************************************************************
      SUBROUTINE PROPAGATOR123(alambda,X1,Y1,Z1,X2,Y2,Z2,d,&
     &CTXX,CTXY,CTXZ,CTYX,CTYY,CTYZ,CTZX,CTZY,CTZZ,&
     &cn1,cn2,cn3)
     
      USE PYGDMPRECISION, ONLY : dp
      implicit none
      
      complex(dp), parameter :: c0=(0._dp,0._dp)
      complex(dp), parameter :: cun=(1._dp,0._dp)
      complex(dp), parameter :: cim=(0._dp,1._dp)
      real(dp), parameter :: Pi=3.141592654_dp
      
      
      complex(dp), intent(in) :: cn1,cn2,cn3
      real(dp), intent(in) :: alambda
      real(dp), intent(in) :: X1,Y1,Z1,X2,Y2,Z2
      real(dp), intent(in) :: d !d: spacing between substrate and top
      
      complex(dp), intent(out) :: CTXX,CTXY,CTXZ,CTYX,CTYY,CTYZ,CTZX,CTZY,CTZZ
      
      
      real(dp) :: ak0
      complex(dp) :: ceps1,ceps2,ceps3
      complex(dp) :: cak1,cak2,cak3,CAK12,CAK22,CAK32
      real(dp) :: ZD,GZ  !d: spacing between substrate and top
!       real(dp) :: TXX1,TYY1,TZZ1,TXY1,TXZ1,TYZ1,TYX1,TZX1,TZY1
!       real(dp) :: TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2,TYX2,TZX2,TZY2
      complex(dp) :: cscreen12,cscreen23
      complex(dp) :: cdelta12,cdelta23
      real(dp) :: XC,YC,ZC
      real(dp) :: T1XX,T1XY,T1XZ,T1YY,T1YZ,T1ZZ
      real(dp) :: T2XX,T2XY,T2XZ,T2YY,T2YZ,T2ZZ
      real(dp) :: T3XX,T3XY,T3XZ,T3YY,T3YZ,T3ZZ
      real(dp) :: RMO
      complex(dp) :: CFEXP
      complex(dp) :: SXX12,SXY12,SXZ12,SYX12,SYY12,SYZ12,SZX12,SZY12,SZZ12
      complex(dp) :: SXX23,SXY23,SXZ23,SYX23,SYY23,SYZ23,SZX23,SZY23,SZZ23
      real(dp) :: X,Y,Z
      
      
      
      
     
!*
!************* Optical indexes ********************************
    
      ceps1=cn1*cn1
      ceps2=cn2*cn2
      ceps3=cn3*cn3
      cscreen12=2._dp/(ceps1+ceps2)
      cscreen23=2._dp/(ceps2+ceps3)
      cdelta12=(ceps1-ceps2)/(ceps1+ceps2)
      cdelta23=(ceps3-ceps2)/(ceps3+ceps2)
!************* Wave vectors ***********************************
      ak0=2._dp*Pi/alambda
      cak1=cn1*ak0
      cak2=cn2*ak0
      cak3=cn3*ak0
!*
!*
!***************** Computation of the propagator ***************
!****************** in each medium *****************************
!*
!***************************************************************************
!****************** First slab (1) ** Z2 < 0 (substrate) *******************
!***************************************************************************
      if (Z2.lt.0.) then
      XC=X2-X1
      YC=Y2-Y1
      ZC=Z2-Z1
!C--------------- Free space propagator -----------------------------
!*
      T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      CAK12=cak1**2
      RMO=SQRT(XC**2+YC**2+ZC**2)
      CFEXP=EXP(CIM*cak1*RMO)
      CTXX=CFEXP*(T1XX-CIM*cak1*T2XX-CAK12*T3XX)*cscreen12
      CTYY=CFEXP*(T1YY-CIM*cak1*T2YY-CAK12*T3YY)*cscreen12
      CTZZ=CFEXP*(T1ZZ-CIM*cak1*T2ZZ-CAK12*T3ZZ)*cscreen12
      CTXY=CFEXP*(T1XY-CIM*cak1*T2XY-CAK12*T3XY)*cscreen12
      CTXZ=CFEXP*(T1XZ-CIM*cak1*T2XZ-CAK12*T3XZ)*cscreen12
      CTYZ=CFEXP*(T1YZ-CIM*cak1*T2YZ-CAK12*T3YZ)*cscreen12
      CTYX=CTXY
      CTZX=CTXZ
      CTZY=CTYZ
      ELSE
!***************************************************************************
!************** Second slab (2) ** D >= Z2 >= 0 (environment) ****************
!***************************************************************************
      if (Z2.ge.0.and.Z2.le.d) then
      XC=X2-X1
      YC=Y2-Y1
      ZC=Z2-Z1
      ZD=Z2+Z1
      GZ=ZD-2._dp*d
!C--------------- Free space propagator -----------------------------
!*
      T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      CAK22=cak2**2
      RMO=SQRT(XC**2+YC**2+ZC**2)
      CFEXP=EXP(CIM*cak2*RMO)
      CTXX=CFEXP*(T1XX-CIM*cak2*T2XX-CAK22*T3XX)/ceps2
      CTYY=CFEXP*(T1YY-CIM*cak2*T2YY-CAK22*T3YY)/ceps2
      CTZZ=CFEXP*(T1ZZ-CIM*cak2*T2ZZ-CAK22*T3ZZ)/ceps2
      CTXY=CFEXP*(T1XY-CIM*cak2*T2XY-CAK22*T3XY)/ceps2
      CTXZ=CFEXP*(T1XZ-CIM*cak2*T2XZ-CAK22*T3XZ)/ceps2
      CTYZ=CFEXP*(T1YZ-CIM*cak2*T2YZ-CAK22*T3YZ)/ceps2
      CTYX=CTXY
      CTZX=CTXZ
      CTZY=CTYZ     
!C-------------------- Surface propagator (1,2) --------------------
      Z=ZD
      X=XC
      Y=YC
      SXX12=cdelta12*(Z*Z+Y*Y-2_dp*X*X)/(X*X+Y*Y+Z*Z)**2.5
      SYY12=cdelta12*(Z*Z+X*X-2_dp*Y*Y)/(X*X+Y*Y+Z*Z)**2.5
      SZZ12=cdelta12*(2_dp*Z*Z-X*X-Y*Y)/(X*X+Y*Y+Z*Z)**2.5
      SXY12=cdelta12*(-3_dp*X*Y)/(X*X+Y*Y+Z*Z)**2.5
      SXZ12=cdelta12*(3_dp*X*Z)/(X*X+Y*Y+Z*Z)**2.5
      SYZ12=cdelta12*(3_dp*Y*Z)/(X*X+Y*Y+Z*Z)**2.5
      SYX12=SXY12
      SZX12=-SXZ12
      SZY12=-SYZ12
!*******************************************************************
!*
!C-------------------- Surface propagator (2,3) --------------------
      Z=GZ
      SXX23=cdelta23*(Z*Z+Y*Y-2_dp*X*X)/(X*X+Y*Y+Z*Z)**2.5
      SYY23=cdelta23*(Z*Z+X*X-2_dp*Y*Y)/(X*X+Y*Y+Z*Z)**2.5
      SZZ23=cdelta23*(2_dp*Z*Z-X*X-Y*Y)/(X*X+Y*Y+Z*Z)**2.5
      SXY23=cdelta23*(-3_dp*X*Y)/(X*X+Y*Y+Z*Z)**2.5
      SXZ23=cdelta23*(3_dp*X*Z)/(X*X+Y*Y+Z*Z)**2.5
      SYZ23=cdelta23*(3_dp*Y*Z)/(X*X+Y*Y+Z*Z)**2.5
      SYX23=SXY23
      SZX23=-SXZ23
      SZY23=-SYZ23
!*******************************************************************

      CTXX=CTXX+SXX12+SXX23
      CTYY=CTYY+SYY12+SYY23
      CTZZ=CTZZ+SZZ12+SZZ23
      CTXY=CTXY+SXY12+SXY23
      CTXZ=CTXZ+SXZ12+SXZ23
      CTYZ=CTYZ+SYZ12+SYZ23
      CTYX=CTYX+SYX12+SYX23
      CTZX=CTZX+SZX12+SZX23
      CTZY=CTZY+SZY12+SZY23
!***************************************************************************
!***************** Third slab (3) ** Z2 > D (top cladding) *****************
!***************************************************************************
      ELSE
      XC=X2-X1
      YC=Y2-Y1
      ZC=Z2-Z1
      ZD=Z2+Z1
!C--------------- Free space propagator -----------------------------
!*
      T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C-------------------------------------------------------------------
      T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
      T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
      T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
      T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
!C------------------------------------------------------------------
      CAK32=cak3**2
      RMO=SQRT(XC**2+YC**2+ZC**2)
      CFEXP=EXP(CIM*cak3*RMO)
      CTXX=CFEXP*(T1XX-CIM*cak3*T2XX-CAK32*T3XX)*cscreen23
      CTYY=CFEXP*(T1YY-CIM*cak3*T2YY-CAK32*T3YY)*cscreen23
      CTZZ=CFEXP*(T1ZZ-CIM*cak3*T2ZZ-CAK32*T3ZZ)*cscreen23
      CTXY=CFEXP*(T1XY-CIM*cak3*T2XY-CAK32*T3XY)*cscreen23
      CTXZ=CFEXP*(T1XZ-CIM*cak3*T2XZ-CAK32*T3XZ)*cscreen23
      CTYZ=CFEXP*(T1YZ-CIM*cak3*T2YZ-CAK32*T3YZ)*cscreen23
      CTYX=CTXY
      CTZX=CTXZ
      CTZY=CTYZ
      ENDIF
      ENDIF
      END

      
      
! !C******************************************************************
! !C****************** propagator_123.f ******************************
! !C******************************************************************
! !C****************** Author: Ch. Girard ****************************
! !C****************** Date: 9 April 1997 ****************************
! !C****************** Based on the retarded image theory ************
! !C******************************************************************
!       SUBROUTINE PROPAGATOR123(alambda,X1,Y1,Z1,X2,Y2,Z2,d,&
!      &TXX1,TYY1,TZZ1,&
!      &TXY1,TXZ1,TYZ1,TYX1,TZX1,TZY1,&
!      &TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2,TYX2,TZX2,TZY2,&
!      &cn1,cn2,cn3)
!      
!       USE PYGDMPRECISION, ONLY : dp
!       implicit none
!       
!       complex(dp), parameter :: c0=(0._dp,0._dp)
!       complex(dp), parameter :: cun=(1._dp,0._dp)
!       complex(dp), parameter :: cim=(0._dp,1._dp)
!       real(dp), parameter :: Pi=3.141592654_dp
!       
!       
!       complex(dp) :: cn1,cn2,cn3
!       complex(dp) :: ceps1,ceps2,ceps3
!       real(dp) :: alambda,ak0
!       complex(dp) :: cak1,cak2,cak3,CAK12,CAK22,CAK32
!       real(dp) :: X1,Y1,Z1,X2,Y2,Z2,d,ZD,GZ  !d: spacing between substrate and top
!       real(dp) :: TXX1,TYY1,TZZ1
!       real(dp) :: TXY1,TXZ1,TYZ1,TYX1,TZX1,TZY1
!       real(dp) :: TXX2,TYY2,TZZ2,TXY2,TXZ2,TYZ2,TYX2,TZX2,TZY2
!       complex(dp) :: cscreen12,cscreen23
!       complex(dp) :: cdelta12,cdelta23
!       real(dp) :: XC,YC,ZC
!       real(dp) :: T1XX,T1XY,T1XZ,T1YY,T1YZ,T1ZZ
!       real(dp) :: T2XX,T2XY,T2XZ,T2YY,T2YZ,T2ZZ
!       real(dp) :: T3XX,T3XY,T3XZ,T3YY,T3YZ,T3ZZ
!       real(dp) :: RMO
!       complex(dp) :: CFEXP
!       complex(dp) :: CTXX,CTXY,CTXZ,CTYX,CTYY,CTYZ,CTZX,CTZY,CTZZ
!       complex(dp) :: SXX12,SXY12,SXZ12,SYX12,SYY12,SYZ12,SZX12,SZY12,SZZ12
!       complex(dp) :: SXX23,SXY23,SXZ23,SYX23,SYY23,SYZ23,SZX23,SZY23,SZZ23
!       real(dp) :: X,Y,Z
!       
!       
!       
!       
!      
! !*
! !************* Optical indexes ********************************
!     
!       ceps1=cn1*cn1
!       ceps2=cn2*cn2
!       ceps3=cn3*cn3
!       cscreen12=2._dp/(ceps1+ceps2)
!       cscreen23=2._dp/(ceps2+ceps3)
!       cdelta12=(ceps1-ceps2)/(ceps1+ceps2)
!       cdelta23=(ceps3-ceps2)/(ceps3+ceps2)
! !************* Wave vectors ***********************************
!       ak0=2._dp*Pi/alambda
!       cak1=cn1*ak0
!       cak2=cn2*ak0
!       cak3=cn3*ak0
! !*
! !*
! !***************** Computation of the propagator ***************
! !****************** in each medium *****************************
! !*
! !****************** First slab (1) ** Z2 < 0 *******************
!       if (Z2.lt.0.) then
!       XC=X2-X1
!       YC=Y2-Y1
!       ZC=Z2-Z1
! !C--------------- Free space propagator -----------------------------
! !*
!       T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       CAK12=cak1**2
!       RMO=SQRT(XC**2+YC**2+ZC**2)
!       CFEXP=EXP(CIM*cak1*RMO)
!       CTXX=CFEXP*(T1XX-CIM*cak1*T2XX-CAK12*T3XX)*cscreen12
!       CTYY=CFEXP*(T1YY-CIM*cak1*T2YY-CAK12*T3YY)*cscreen12
!       CTZZ=CFEXP*(T1ZZ-CIM*cak1*T2ZZ-CAK12*T3ZZ)*cscreen12
!       CTXY=CFEXP*(T1XY-CIM*cak1*T2XY-CAK12*T3XY)*cscreen12
!       CTXZ=CFEXP*(T1XZ-CIM*cak1*T2XZ-CAK12*T3XZ)*cscreen12
!       CTYZ=CFEXP*(T1YZ-CIM*cak1*T2YZ-CAK12*T3YZ)*cscreen12
!       CTYX=CTXY
!       CTZX=CTXZ
!       CTZY=CTYZ
!       TXX1=REAL(CTXX)
!       TYY1=REAL(CTYY)
!       TZZ1=REAL(CTZZ)
!       TXY1=REAL(CTXY)
!       TXZ1=REAL(CTXZ)
!       TYZ1=REAL(CTYZ)
!       TYX1=REAL(CTYX)
!       TZX1=REAL(CTZX)
!       TZY1=REAL(CTZY)
!       TXX2=AIMAG(CTXX)
!       TYY2=AIMAG(CTYY)
!       TZZ2=AIMAG(CTZZ)
!       TXY2=AIMAG(CTXY)
!       TXZ2=AIMAG(CTXZ)
!       TYZ2=AIMAG(CTYZ)
!       TYX2=AIMAG(CTYX)
!       TZX2=AIMAG(CTZX)
!       TZY2=AIMAG(CTZY)
!       ELSE
!       if (Z2.ge.0.and.Z2.le.d) then
!       XC=X2-X1
!       YC=Y2-Y1
!       ZC=Z2-Z1
!       ZD=Z2+Z1
!       GZ=ZD-2._dp*d
! !C--------------- Free space propagator -----------------------------
! !*
!       T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       CAK22=cak2**2
!       RMO=SQRT(XC**2+YC**2+ZC**2)
!       CFEXP=EXP(CIM*cak2*RMO)
!       CTXX=CFEXP*(T1XX-CIM*cak2*T2XX-CAK22*T3XX)/ceps2
!       CTYY=CFEXP*(T1YY-CIM*cak2*T2YY-CAK22*T3YY)/ceps2
!       CTZZ=CFEXP*(T1ZZ-CIM*cak2*T2ZZ-CAK22*T3ZZ)/ceps2
!       CTXY=CFEXP*(T1XY-CIM*cak2*T2XY-CAK22*T3XY)/ceps2
!       CTXZ=CFEXP*(T1XZ-CIM*cak2*T2XZ-CAK22*T3XZ)/ceps2
!       CTYZ=CFEXP*(T1YZ-CIM*cak2*T2YZ-CAK22*T3YZ)/ceps2
!       CTYX=CTXY
!       CTZX=CTXZ
!       CTZY=CTYZ     
! !C-------------------- Surface propagator (1,2) --------------------
!       Z=ZD
!       X=XC
!       Y=YC
!       SXX12=cdelta12*(Z*Z+Y*Y-2_dp*X*X)/(X*X+Y*Y+Z*Z)**2.5
!       SYY12=cdelta12*(Z*Z+X*X-2_dp*Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SZZ12=cdelta12*(2_dp*Z*Z-X*X-Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SXY12=cdelta12*(-3_dp*X*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SXZ12=cdelta12*(3_dp*X*Z)/(X*X+Y*Y+Z*Z)**2.5
!       SYZ12=cdelta12*(3_dp*Y*Z)/(X*X+Y*Y+Z*Z)**2.5
!       SYX12=SXY12
!       SZX12=-SXZ12
!       SZY12=-SYZ12
! !*******************************************************************
! !*
! !C-------------------- Surface propagator (2,3) --------------------
!       Z=GZ
!       SXX23=cdelta23*(Z*Z+Y*Y-2_dp*X*X)/(X*X+Y*Y+Z*Z)**2.5
!       SYY23=cdelta23*(Z*Z+X*X-2_dp*Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SZZ23=cdelta23*(2_dp*Z*Z-X*X-Y*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SXY23=cdelta23*(-3_dp*X*Y)/(X*X+Y*Y+Z*Z)**2.5
!       SXZ23=cdelta23*(3_dp*X*Z)/(X*X+Y*Y+Z*Z)**2.5
!       SYZ23=cdelta23*(3_dp*Y*Z)/(X*X+Y*Y+Z*Z)**2.5
!       SYX23=SXY23
!       SZX23=-SXZ23
!       SZY23=-SYZ23
! !*******************************************************************
! 
!       TXX1=REAL(CTXX+SXX12+SXX23)
!       TYY1=REAL(CTYY+SYY12+SYY23)
!       TZZ1=REAL(CTZZ+SZZ12+SZZ23)
!       TXY1=REAL(CTXY+SXY12+SXY23)
!       TXZ1=REAL(CTXZ+SXZ12+SXZ23)
!       TYZ1=REAL(CTYZ+SYZ12+SYZ23)
!       TYX1=REAL(CTYX+SYX12+SYX23)
!       TZX1=REAL(CTZX+SZX12+SZX23)
!       TZY1=REAL(CTZY+SZY12+SZY23)
!       TXX2=AIMAG(CTXX+SXX12+SXX23)
!       TYY2=AIMAG(CTYY+SYY12+SYY23)
!       TZZ2=AIMAG(CTZZ+SZZ12+SZZ23)
!       TXY2=AIMAG(CTXY+SXY12+SXY23)
!       TXZ2=AIMAG(CTXZ+SXZ12+SXZ23)
!       TYZ2=AIMAG(CTYZ+SYZ12+SYZ23)
!       TYX2=AIMAG(CTYX+SYX12+SYX23)
!       TZX2=AIMAG(CTZX+SZX12+SZX23)
!       TZY2=AIMAG(CTZY+SZY12+SZY23)
!       ELSE
!       XC=X2-X1
!       YC=Y2-Y1
!       ZC=Z2-Z1
!       ZD=Z2+Z1
! !C--------------- Free space propagator -----------------------------
! !*
!       T1XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XX=(2_dp*XC*XC-YC*YC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XX=-(YC*YC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XY=3_dp*XC*YC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XY=XC*YC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2XZ=3_dp*XC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3XZ=XC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YY=(2_dp*YC*YC-XC*XC-ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YY=-(XC*XC+ZC*ZC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C-------------------------------------------------------------------
!       T1YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2YZ=3_dp*YC*ZC/(XC*XC+YC*YC+ZC*ZC)**2
!       T3YZ=YC*ZC/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       T1ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2.5
!       T2ZZ=(2_dp*ZC*ZC-XC*XC-YC*YC)/(XC*XC+YC*YC+ZC*ZC)**2
!       T3ZZ=-(XC*XC+YC*YC)/(XC*XC+YC*YC+ZC*ZC)**1.5
! !C------------------------------------------------------------------
!       CAK32=cak3**2
!       RMO=SQRT(XC**2+YC**2+ZC**2)
!       CFEXP=EXP(CIM*cak3*RMO)
!       CTXX=CFEXP*(T1XX-CIM*cak3*T2XX-CAK32*T3XX)*cscreen23
!       CTYY=CFEXP*(T1YY-CIM*cak3*T2YY-CAK32*T3YY)*cscreen23
!       CTZZ=CFEXP*(T1ZZ-CIM*cak3*T2ZZ-CAK32*T3ZZ)*cscreen23
!       CTXY=CFEXP*(T1XY-CIM*cak3*T2XY-CAK32*T3XY)*cscreen23
!       CTXZ=CFEXP*(T1XZ-CIM*cak3*T2XZ-CAK32*T3XZ)*cscreen23
!       CTYZ=CFEXP*(T1YZ-CIM*cak3*T2YZ-CAK32*T3YZ)*cscreen23
!       CTYX=CTXY
!       CTZX=CTXZ
!       CTZY=CTYZ
!       TXX1=REAL(CTXX)
!       TYY1=REAL(CTYY)
!       TZZ1=REAL(CTZZ)
!       TXY1=REAL(CTXY)
!       TXZ1=REAL(CTXZ)
!       TYZ1=REAL(CTYZ)
!       TYX1=REAL(CTYX)
!       TZX1=REAL(CTZX)
!       TZY1=REAL(CTZY)
!       TXX2=AIMAG(CTXX)
!       TYY2=AIMAG(CTYY)
!       TZZ2=AIMAG(CTZZ)
!       TXY2=AIMAG(CTXY)
!       TXZ2=AIMAG(CTXZ)
!       TYZ2=AIMAG(CTYZ)
!       TYX2=AIMAG(CTYX)
!       TZX2=AIMAG(CTZX)
!       TZY2=AIMAG(CTZY)
!       ENDIF
!       ENDIF
!       END
