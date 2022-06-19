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
! Date: March 1997
!
! Adapted for python-interface 10/2017 by Peter Wiecha
!
!*********************************************************************
!*********************************************************************
! - Double interfaces located at z = 0 and z = D 
!
! - Warning: Problem on sign of magnetic field 
!*********************************************************************



    subroutine evanescentfield(alambda,thetai,d,x,y,z, &
                                        & cn1,cn2,cn3, &
                                        & cex,cey,cez, &
                                        & cmagx,cmagy,cmagz)
!**************************************************************
    USE PYGDMPRECISION, ONLY : dp
    implicit none
    
    
    complex(dp), parameter :: C0=(0._dp,0._dp)
    complex(dp), parameter :: CUN=(1._dp,0._dp)
    complex(dp), parameter :: CIM=(0._dp,1._dp)
    real(dp), parameter :: Pi=3.141592654_dp
    real(dp), parameter :: e0=1.0_dp
    real(dp), parameter :: b0=1.0_dp
!**************************************************************
    
    real(dp), intent(in) :: alambda,thetai,d,x,y,z
    complex(dp), intent(in) :: cn1, cn2, cn3
    
    complex(dp), intent(out) :: cex(2),cey(2),cez(2)
    complex(dp), intent(out) :: cmagx(2),cmagy(2),cmagz(2)

    
!**************************************************************
    real(dp) :: ak0, theta1
    
    complex(dp) :: ceps1
    complex(dp) :: ceps2
    complex(dp) :: ceps3

    complex(dp) :: ceffec12, ceffec13
    complex(dp) :: cak1x, cak1y, cak1z, cagk1z
    complex(dp) :: cak2x, cak2y, cak2z, cagk2z
    complex(dp) :: cak3x, cak3y, cak3z, cagk3z
    complex(dp) :: cakp1x, cakp1y, cakp1z, cagkp1z
    complex(dp) :: cakp2x, cakp2y, cakp2z, cagkp2z
    
    complex(dp) :: c2p, c2m, c3
    complex(dp) :: cdelta, cdelta1p, cdelta2, cdelta2p, cdelta3
    complex(dp) :: cep1, ce2, cep2, ce3
    complex(dp) :: cdeltatm, cdelta1ptm, cdelta2tm, cdelta2ptm, cdelta3tm
    complex(dp) :: cmagp1, cmag2, cmagp2, cmag3
    complex(dp) :: cphase1, cphase1p, cphase2, cphase2p, cphase3
                           
    
!************* Incident angle in radian ***********************
    theta1=thetai*Pi/180._dp
    
!************* Optical indexes ********************************
    ceps1=cn1*cn1
    ceps2=cn2*cn2
    ceps3=cn3*cn3

!************* Wave vectors ***********************************
    ak0=2.*Pi/alambda
    
    ceffec12 = cun*sqrt(ceps2-ceps1*sin(theta1) &
              & * sin(theta1))
    ceffec13 = cun*sqrt(ceps3-ceps1*sin(theta1) &
              & * sin(theta1))
!*
!*
!************* First medium (1) *******************************
!*
!****** Incident wave *****************************************
    cak1x=-ak0*cn1*sin(theta1)
    cak1y=c0
    cak1z=ak0*cn1*cos(theta1)
    cagk1z=cak1z/ceps1
!****** Reflected wave ****************************************
    cakp1x=-ak0*cn1*sin(theta1)
    cakp1y=c0
    cakp1z=-ak0*cn1*cos(theta1) 
    cagkp1z=cakp1z/ceps1
!************* Second medium (2) ******************************
!*
!******* First wave ******************************************* 
    cak2x=-ak0*cn1*sin(theta1)
    cak2y=c0
    cak2z=cun*ak0*ceffec12
    cagk2z=cak2z/ceps2
!*
!******* Second wave ******************************************
    cakp2x=-ak0*cn1*sin(theta1)
    cakp2y=c0
    cakp2z=-cun*ak0*ceffec12
    cagkp2z=cakp2z/ceps2
!**************************************************************
!*
!************** Third medium **********************************
    cak3x=-ak0*cn1*sin(theta1)
    cak3y=c0
    cak3z=cun*ak0*ceffec13
    cagk3z=cak3z/ceps3
!*
!************** Phase factors *********************************
    c2p=exp(cim*cak2z*d)
    c2m=exp(-cim*cak2z*d)
    c3=exp(cim*cak3z)
!*
!C#############################################################
!*
!*
!*********** Computation of the electric field ****************
!******* modulus in the case of the s-polarized mode **********
!**************************************************************
!*
!*************** Main determinant *****************************
    cdelta = (c2p*(cak2z*cak2z-cak3z*cak2z+cak1z*cak3z - &
              & cak1z*cak2z)-c2m*(cak2z*cak2z+cak3z*cak2z + &
              & cak1z*cak3z+cak1z*cak2z))
!**************************************************************
!*
!*********** Computation of the other determinants ************
!*
!*********************** Ep1 **********************************
    cdelta1p = e0*(c2p*(-cak2z*cak2z+cak3z*cak2z+cak1z*cak3z - &
                & cak1z*cak2z)+c2m*(cak2z*cak2z+cak3z*cak2z - &
                & cak1z*cak3z-cak1z*cak2z))
!*
!*********************** E2 ***********************************
    cdelta2=-e0*c2m*(2.*cak1z*cak3z+2.*cak1z*cak2z)
!*
!*********************** Ep2 **********************************
    cdelta2p=e0*c2p*(-2.*cak1z*cak2z+2.*cak1z*cak3z)
!*
!*********************** E3 ***********************************
    cdelta3=-4.*e0*c2p*c2m*cak1z*cak2z
!*
!**************************************************************
!*
!********** Complex amplitudes inside the three media *********
!*
!*************** (1) Electric field ***************************
    cep1=cdelta1p/cdelta
    ce2=cdelta2/cdelta
    cep2=cdelta2p/cdelta
    ce3=cdelta3/(cdelta*c3)
    
!**************************************************************
!*
!C#############################################################
!*
!*
!*********** Computation of the magnetic field ****************
!******* modulus in the case of the p-polarized mode **********
!**************************************************************
!*
!*************** Main determinant *****************************
    cdeltatm = (c2p*(cagk2z*cagk2z-cagk3z*cagk2z+cagk1z*cagk3z - &
                & cagk1z*cagk2z)-c2m*(cagk2z*cagk2z+cagk3z*cagk2z + &
                & cagk1z*cagk3z+cagk1z*cagk2z))
!**************************************************************
!*
!*********** Computation of the other determinants ************
!*
!*********************** Bp1 **********************************
    cdelta1ptm = e0*(c2p*(-cagk2z*cagk2z+cagk3z*cagk2z + &
                & cagk1z*cagk3z - &
                & cagk1z*cagk2z)+c2m*(cagk2z*cagk2z+cagk3z*cagk2z - &
                & cagk1z*cagk3z-cagk1z*cagk2z))
!*
!*********************** B2 ***********************************
    cdelta2tm=-e0*c2m*(2.*cagk1z*cagk3z+2.*cagk1z*cagk2z)
!*
!*********************** Bp2 **********************************
    cdelta2ptm=e0*c2p*(-2.*cagk1z*cagk2z+2.*cagk1z*cagk3z)
!*
!*********************** B3 ***********************************
    cdelta3tm=-4.*e0*c2p*c2m*cagk1z*cagk2z
!*
!**************************************************************
!*
!********** Complex amplitudes inside the three media *********
!*
!*
!*************** (2) Magnetic field ***************************
    cmagp1 = cdelta1ptm / cdeltatm
    cmag2  = cdelta2tm / cdeltatm
    cmagp2 = cdelta2ptm / cdeltatm
    cmag3  = cdelta3tm / (cdeltatm*c3)

!**************************************************************
!*** Electric and Magnetic fields inside each material slab ***
!*
!****************** First medium (z < 0 )**********************
    if (z.lt.0.) then
        cphase1=exp(cim*(cak1x*x+cak1y*y+cak1z*z))
        cphase1p=exp(cim*(cakp1x*x+cakp1y*y+cakp1z*z))
!***************** Elec ***************************************
        cex(1)=c0
        cey(1)=e0*cphase1+cep1*cphase1p
        cez(1)=c0
!*-----------------Mag ----------------------------------------
        cmagx(1) = -(e0*cphase1*cak1z/ak0) &
                 & - (cep1*cphase1p*cakp1z/ak0)
        cmagy(1)=c0
        cmagz(1)=(e0*cphase1*cak1x/ak0)+(cep1*cphase1p*cakp1x/ak0)
!*-------------------------------------------------------------
!*
!*----------------- p-polarized mode --------------------------
        cmagx(2)=-c0
        cmagy(2)=-(b0*cphase1+cmagp1*cphase1p)
        cmagz(2)=-c0
        cex(2) = (cn1)*(-b0*cphase1*cak1z/(ceps1*ak0) &
                & -cmagp1*cphase1p*cakp1z/(ceps1*ak0))
        cey(2)=c0
        cez(2) = (cn1)*(b0*cphase1*cak1x/(ceps1*ak0) &
                & +cmagp1*cphase1p*cakp1x/(ceps1*ak0))
!*-------------------------------------------------------------
    else
        if (z.ge.0.and.z.le.d) then
            cphase2=exp(cim*(cak2x*x+cak2y*y+cak2z*z))
            cphase2p=exp(cim*(cakp2x*x+cakp2y*y+cakp2z*z))
            cex(1)=c0
            cey(1)=ce2*cphase2+cep2*cphase2p
            cez(1)=c0
!*
!*-------------------------------------------------------------
            cmagx(1) =  -(ce2*cphase2*cak2z/ak0) &
                      & -(cep2*cphase2p*cakp2z/ak0)
            cmagy(1)=c0
            cmagz(1)=(ce2*cphase2*cak2x/ak0)+(cep2*cphase2p*cakp2x/ak0)
!*-------------------------------------------------------------
!*----------------- p-polarized mode --------------------------
            cmagx(2)=-c0
            cmagy(2)=-(cmag2*cphase2+cmagp2*cphase2p)
            cmagz(2)=-c0
            cex(2) = (cn1)*(-cmag2*cphase2*cak2z/(ceps2*ak0) &
                    & -cmagp2*cphase2p*cakp2z/(ceps2*ak0))
            cey(2)=c0
            cez(2) = (cn1)*(cmag2*cphase2*cak2x/(ceps2*ak0) &
                    & +cmagp2*cphase2p*cakp2x/(ceps2*ak0))
        
        else
!*--------------------------------------------------------------
            cphase3=exp(cim*(cak3x*x+cak3y*y+cak3z*z))
            cex(1)=c0
            cey(1)=ce3*cphase3
            cez(1)=c0
!*-------------------------------------------------------------
            cmagx(1)=-(ce3*cphase3*cak3z/ak0)
            cmagy(1)=c0
            cmagz(1)=(ce3*cphase3*cak3x/ak0)
!*----------------- p-polarized mode --------------------------
            cmagx(2)=-c0
            cmagy(2)=-cmag3*cphase3
            cmagz(2)=-c0
            cex(2)=(cn1)*(-cmag3*cphase3*cak3z/(ceps3*ak0))
            cey(2)=c0
            cez(2)=(1./cn1)*(cmag3*cphase3*cak3x/(ceps3*ak0))
        endif
    endif
    end SUBROUTINE evanescentfield
      


      
