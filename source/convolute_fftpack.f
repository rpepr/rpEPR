! An f77 subroutine to convolute a stick spectrum with a derivative
! gaussian line shape using fftpack5 routines
      subroutine convolute(spec,n,lw,sw)
      implicit none
      !Argument variables
      real spec(10000),lw,sw
      integer n

      !Local variables
      real wsave(33000),work(65536),r, pi
      complex scr1(32768),myi
      integer lensave,ierr,lenspec,lenwork,i,kb
      external cfft1i
      external cfft1f
      external cfft1b

      pi=3.14159265359
      kb = 1024
      myi = cmplx(0,1)
      lensave = 33000
      lenwork = 655536

! zero fill the spectrum to the next nearest power of 2
      do i = 1,5
        if(kb .gt. n) goto 10
        kb = kb * 2
      end do
 10   lenspec = kb
      do i=1,n
        scr1(i+(kb-n)/2) = cmplx(spec(i),0)
      end do

! Set up and perform the fft
      call cfft1i(kb,wsave,lensave,ierr)
      call cfft1f(kb,1,scr1,lenspec,wsave,lensave,work,lenwork,ierr)

! Multiply the FT by the convolution function
      r = 0.
      do i = 1, kb
        scr1(i) = scr1(i)*exp(-.5*(pi*lw*r/sw)**2) ! FT of a gaussian
        scr1(i) = scr1(i)*2.*pi*r*myi/n ! time domain first derivative
        r = r + 1.
      end do

! Perform the reverse FFT
      call cfft1b(kb,1,scr1,lenspec,wsave,lensave,work,lenwork,ierr)
      do i = 1,n
        spec(i) = n*real(scr1(i+(kb-n)/2))
      end do
      end subroutine convolute
