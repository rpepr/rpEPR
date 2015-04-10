      module gaussian_line
      implicit none
      contains
!     ***********************************************************
!     Subroutine to calculate a gaussian lineshape
!     ***********************************************************

      subroutine gline(lshape,lw,xstep,lpts)
      implicit none
      real :: lshape(10000)
      real :: lw, xstep, den, X

      integer lpts,i

      den = lw**2.
      do i= 1, 10000
         lshape(i) = 0
      end do

      X = -5*lw
      do i=1,lpts
         lshape(i) = (0.47/lw)*(-1.38*X/den)*exp(((-0.693)*X**2.)/den)
         lshape(i) = lshape(i) * 100.
         X = X + xstep
      end do

      return
      END subroutine gline
      end module gaussian_line
