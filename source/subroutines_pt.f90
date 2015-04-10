       module setup_routines
       implicit none
       contains
!     **********************************************************
!     Subrotine to do Mt*N*O on 3x3 matrices, result is returned
!     in the middle matrix.
!     *********************************************************

      subroutine multran(A,B,C)
      real A(3,3), B(3,3), C(3,3), S(3,3)
      integer i,j,k,l

      do i=1,3
         do j=1,3
            S(i,j)=0
         end do
      end do

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  S(i,j)=S(i,j)+A(k,i)*B(k,l)*C(l,j)
               end do
            end do
         end do
      end do

      do i=1,3
         do j=1,3
            B(i,j)=S(i,j)
         end do
      end do

      return
      END subroutine multran

!     ***********************************************************
!     Subroutine to rotate a tensor using a rotation matrix
!     ***********************************************************

      subroutine rotate(D,R)
!     D is the 3x3 tensor to be rotated, and R is the rotation
!     matrix, S is the rotated matrix

      real D(3,3),R(3,3),S(3,3)
      integer i,j,k,l

      do j=1,3
         do i=1,3
            S(i,j)=0.
         end do
      end do

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  S(i,j)=S(i,j)+R(i,k)*D(k,l)*R(j,l)
               end do
            end do
         end do
      end do

      do j=1,3
         do i=1,3
            D(i,j)=S(i,j)
         end do
      end do

      return
      END subroutine rotate
      end module setup_routines


      module spectrum_calculate
      implicit none
      contains
!     ***********************************************************
!     Function to calculate the weighted average of two resonance
!     estimates.
!     ***********************************************************

      real function wav_reson(field1,field2,res_est1,res_est2)
      real field1,field2,res_est1,res_est2

      wav_reson = (res_est1*field2 - res_est2*field1)/(field2 -         &
     &             field1 + res_est1 - res_est2)

!      res_est = (res_est1+res_est2)/2.
!      do i = 1,20
!         wgt1 = abs(res_est-field2)
!         wgt2 = abs(res_est-field1)
!         wav_reson = (res_est1*wgt1+res_est2*wgt2)/(wgt1+wgt2)
!         if (abs(res_est-wav_reson).lt.tol) return
!         res_est = wav_reson
!      end do
!      write(*,*) "Error calculating weighted average resonance."
      return
      END function wav_reson

!     **********************************************************
!     Subroutine to add a pre-calculated line to a spectrum
!     **********************************************************

      subroutine lineadd(H,HR,lshape,TP,xstep,hlow,lpts)

      real, dimension(10000) :: H, lshape

      real :: HR, xstep, hlow

      real TP

      integer respt, midpt, step, lpts, clip,i

      midpt = lpts/2
      step = 0

      respt = int((HR-hlow)/xstep)
      clip = 0
      if (midpt .ge. respt) clip=midpt-respt

      do i=respt-midpt+clip+1,respt+midpt
         step = step +1
         H(i) = H(i) + lshape(step+clip)*TP
      end do
      return
      END subroutine lineadd

!     **********************************************************
!     new inproved subroutine to calculate the energy matrix
!     **********************************************************

      subroutine energy(F,Jex,d0,rd1,idp1,idm1,rd2,idp2,idm2,           &
     & rd3,idp3,idm3,rd4,idp4,idm4,nucqnum,nnuc,ne,ar,ai,               &
     & Atens,gnlab,glab1,glab2)

     use atom_types, only: maxatoms

      real ar(4,4), ai(4,4),bohrmag,                &
     & qnum,d0,rd1,idp1,idm1,rd2,idp2,idm2,                             &
     & rd3,idp3,idm3,rd4,idp4,idm4,nucmag,                              &
     & r1z0, r2z0, r1z1, r2z1, i1zp1, i2zp1, i1zm1, i2zm1,              &
     & a0, ra1, iap1, iam1,glab1(3,3),glab2(3,3),                       &
     & Jex,Atens(maxatoms,3,3),gnlab(maxatoms,3,3)


      real F,nucqnum(maxatoms)

      integer ne(maxatoms),nnuc,k,i,j

      bohrmag=9.2740154e-21
      nucmag=5.0507866e-24

      do i=1,4
         do j=1,4
            ar(i,j) = 0.
            ai(i,j) = 0.
         end do
      end do

!      print*,'d0',d0,'rd1',rd1,'idp1',idp1,'idm1',
!     .       idm1,'rd2',rd2,
!     .       'idp2',idp2,'idm2',idm2,'rd3',rd3,'idp3',
!     .       idp3,'idm3',idm3,
!     .       'rd4',rd4,'idp4',idp4,'idm4',idm4

!     Define Zeeman terms

!      print*,'glab1'
!      print*,glab1(1,1),glab1(1,2),glab1(1,3)
!      print*,glab1(2,1),glab1(2,2),glab1(2,3)
!      print*,glab1(3,1),glab1(3,2),glab1(3,3)

!      print*,'glab2'
!      print*,glab2(1,1),glab2(1,2),glab2(1,3)
!      print*,glab2(2,1),glab2(2,2),glab2(2,3)
!      print*,glab2(3,1),glab2(3,2),glab2(3,3)

      r1z0 = bohrmag*F*glab1(3,3)
      r2z0 = bohrmag*F*glab2(3,3)
      r1z1 = bohrmag*F*glab1(3,1)/2.
      r2z1 = bohrmag*F*glab2(3,1)/2.
      i1zp1 = -bohrmag*F*glab1(3,2)/2.
      i1zm1 = -i1zp1
      i2zp1 = -bohrmag*F*glab2(3,2)/2.
      i2zm1 = -i2zp1

!      print*, F

!      print*, 'rz0',r1z0,'rz1',r1z1,'izp1',i1zp1,'izm1',i1zm1
!      print*, 'rz0',r2z0,'rz1',r2z1,'izp1',i2zp1,'izm1',i2zm1

!     Create energy matrix.

      ar(1,1) = (r1z0 + r2z0)/2. + Jex/4. + d0/4.
      ar(1,2) = r2z1 + rd2/2.
      ai(1,2) = i2zp1 + idp2/2.
      ar(1,3) = r1z1 + rd1/2.
      ai(1,3) = i1zp1 + idp1/2.
      ar(1,4) = rd4
      ai(1,4) = idp4
      ar(2,1) = r2z1 + rd2/2.
      ai(2,1) = i2zm1 + idm2/2.
      ar(2,2) = (r1z0 - r2z0)/2. - Jex/4. - d0/4.
      ar(2,3) = Jex/2. + rd3
      ai(2,3) = idp3
      ar(2,4) = r1z1 - rd1/2.
      ai(2,4) = i1zp1 - idp1/2.
      ar(3,1) = r1z1 + rd1/2.
      ai(3,1) = i1zm1 + idm1/2.
      ar(3,2) = Jex/2. + rd3
      ai(3,2) = idm3
      ar(3,3) = -(r1z0 - r2z0)/2. - Jex/4. - d0/4.
      ar(3,4) = r2z1 - rd2/2.
      ai(3,4) = i2zp1 - idp2/2.
      ar(4,1) = rd4
      ai(4,1) = idm4
      ar(4,2) = r1z1 - rd1/2.
      ai(4,2) = i1zm1 - idm1/2.
      ar(4,3) = r2z1 - rd2/2.
      ai(4,3) = i2zm1 - idm2/2.
      ar(4,4) = (-r1z0 - r2z0)/2. + Jex/4. + d0/4.

!     Add the nuclear hyperfine terms

      do k=1,nnuc
         a0 = Atens(k,3,3)
         ra1 = Atens(k,1,3)/2.
         iap1 = -Atens(k,2,3)/2.
         iam1 = -iap1
         qnum = nucqnum(k)
         if(ne(k) .eq. 1) then
            ar(1,1) = ar(1,1) + qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(1,3) = ar(1,3) + qnum*ra1
            ai(1,3) = ai(1,3) + qnum*iap1
            ar(2,2) = ar(2,2) + qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(2,4) = ar(2,4) + qnum*ra1
            ai(2,4) = ai(2,4) + qnum*iap1
            ar(3,1) = ar(3,1) + qnum*ra1
            ai(3,1) = ai(3,1) + qnum*iam1
            ar(3,3) = ar(3,3) - qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(4,2) = ar(4,2) + qnum*ra1
            ai(4,2) = ai(4,2) + qnum*iam1
            ar(4,4) = ar(4,4) - qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
         else
            ar(1,1) = ar(1,1) + qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(1,2) = ar(1,2) + qnum*ra1
            ai(1,2) = ai(1,2) + qnum*iap1
            ar(2,1) = ar(2,1) + qnum*ra1
            ai(2,1) = ai(2,1) + qnum*iam1
            ar(2,2) = ar(2,2) - qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(3,4) = ar(3,4) + qnum*ra1
            ai(3,4) = ai(3,4) + qnum*iap1
            ar(3,3) = ar(3,3) + qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
            ar(4,3) = ar(4,3) + qnum*ra1
            ai(4,3) = ai(4,3) + qnum*iam1
            ar(4,4) = ar(4,4) - qnum*a0/2.                              &
     &      + qnum*nucmag*F*gnlab(k,3,3)
         end if
      end do

!      do i=1,4
!         print*, (ar(i,j),j=1,4)
!         print*, (ai(i,j),j=1,4)
!      end do

      return
      END subroutine energy

      end module spectrum_calculate
