! An f77 subroutine to calculate the correlation of 2 spectra using fftpack5 routines
      function corrf(spec1,spec2,n)
      implicit none
      !Argument variables
      real spec1(10000),spec2(10000)
      real corrf(10000)
      integer n
      !Local variables
      real mycorr(10000),wsave(10008),work(n),scr1(10000),scr2(10000)
      integer lensave,ierr,lenspec,lenwork,i
      external rfft1i
      external rfft1f
      external rfft1b

      lensave = 10008
      lenspec = 10000
      lenwork = n

      do i=1,n
        scr1(i) = spec1(i)
        scr2(i) = spec2(i)
      end do

! Perform the transforms
      call rfft1i(n,wsave,lensave,ierr)
      call rfft1f(n,1,scr1,lenspec,wsave,lensave,work,lenwork,ierr)
      call rfft1f(n,1,scr2,lenspec,wsave,lensave,work,lenwork,ierr)

! Multiply the FT's together
      do i=1,n,2
        mycorr(i)=scr1(i)*scr2(i)
        mycorr(i+1)=scr1(i+1)*-scr2(i+1)
      end do

! Perform the reverse transform
      call rfft1b(n,1,mycorr,lenspec,wsave,lensave,work,lenwork,ierr)
      corrf = mycorr
      end function
