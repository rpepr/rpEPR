
!     ************************************************************
!     This subroutine calculates the EPR spectrum of a
!     one electron system.
!     ***********************************************************
      module radical_simulation
      contains

      subroutine radsim(spec,H,field,npts)

      use spectrum_parameters
      use atom_types
      use init_atoms
      use gaussian_line
      use physical_constants
      use spectrum_parameter_handlers

      implicit none

      interface
      subroutine convolute(spec,n,lw,sw)
      implicit none
      !Argument variables
      real spec(10000),lw,sw
      integer n
      end subroutine convolute
      end interface

      type(isotope) :: myatom(natoms)
      type(spec_parms), intent(in) :: spec
      type(spec_parms) :: localparms

      real, PARAMETER :: cbeta=9.274015431e-21
      real, parameter :: sq2=1.414213562373
      integer, parameter :: maxres = 10000 ! The size of the resonance storage array

      CHARACTER*5 nucleus(maxatoms)

      real, dimension(10000) :: lshape, H, field, stick
      real :: gps, lwstick, stickrange,specwidth, gmin,gmax,Amax,spinmax,Amin,anmax
      real, allocatable :: nucqnum(:,:), res(:)

      real :: vnu,LW,LWx,LWy,LWz,lwmin,glow,ghigh,Atens(maxatoms,3), &
     &origlow,orighigh,alpha(maxatoms),beta(maxatoms),gam(maxatoms), &
     &xstep,newhigh,newlow,lwmax,ymax,ymin,yshift, &
     &yave,nucspin(maxatoms),hnu,P,geff,dcn,dcm,dcl,S1,S2, &
     &S3,newres

      real, allocatable :: cop(:),sip(:),co2p(:),si2p(:),cot(:), &
     &sit(:),co2t(:),si2t(:)

      real, dimension(maxatoms) :: sial,sibe,siga,coal,cobe,coga,aeff


      INTEGER :: I,R,lpts,npad,nnuc,nucmult,LP,LT, &
     & num,numpts,A,B,npts,j,nres

      integer, allocatable :: istick(:)

      logical :: eulerflag, lwiso

      real :: g(3)

     call fill_atoms(myatom)

     ! Calculate A values for beta nuclei, and equivalent nuclei
      localparms = spec
      call parms_to_A(localparms)

! Copy the spectral parameters into local variables
      nnuc =  localparms%numnuclei
      nucmult=1

      DO i=1,nnuc ! calculate the nuclear multiplicity, and allocate storage for the quantum numbers
        nucleus(i) = localparms%nucname(i)
        nucspin(i) = myatom(localparms%what_isotope(i))%isospin
        nucmult=nucmult*(2*nucspin(i)+1)
      END DO
      allocate(nucqnum(nnuc,nucmult))

      origlow = localparms%origlow
      orighigh = localparms%orighigh

      vnu = localparms%freq
      vnu=1e9*vnu
      hnu=vnu*plank

      LT = localparms%lt
      LP = localparms%lp

      g=0.

      g(1) = localparms%g1(1)
      g(2) = localparms%g1(2)
      g(3) = localparms%g1(3)

      Atens=0.
      eulerflag = .false.
      Do num=1,nnuc
         Atens(num,1) = localparms%Atens(num,1)
         Atens(num,2) = localparms%Atens(num,2)
         Atens(num,3) = localparms%Atens(num,3)
         alpha(num) = localparms%alpha(num)
         beta(num) = localparms%beta(num)
         gam(num) = localparms%gama(num)
         IF (alpha(num) /= 0. .or. beta(num) /= 0. .or. gam(num) /= 0.)eulerflag = .true.
      END DO

      DO num=1,nnuc ! Calculate and store trig functions for euler angles
        alpha(num)=(alpha(num)/180.)*Pi
        beta(num)=(beta(num)/180.)*Pi
        gam(num)=(gam(num)/180.)*Pi
        sial(num) = sin(alpha(num))
        sibe(num) = sin(beta(num))
        siga(num) = sin(gam(num))
        coal(num) = cos(alpha(num))
        cobe(num) = cos(beta(num))
        coga(num) = cos(gam(num))
      END DO

      LWx = localparms%lw(1,1)
      LWy = localparms%lw(1,2)
      LWz = localparms%lw(1,3)
      lwiso=.false.
      IF (LWx.eq.LWy.and.LWy.eq.LWz) lwiso=.true.

!     ***********************************************************
!     The next section calculates and stores values that are independant
!     of the crystal orientation.
!     ***********************************************************

!     a new field range.  We calculate the spectrum over this expanded
!     field range, but only return the portion corresponding to the
!     original field range to avoid artifacts at the edge of the spectrum
!     due to skipping transtions outside the field range.

      lwmax=LWx
      IF (LWy.gt.lwmax) lwmax=LWy
      IF (LWz.gt.lwmax) lwmax=LWz


      newlow=origlow-5*lwmax
      newhigh=orighigh+5*lwmax
      gps=newhigh-newlow
      stickrange = orighigh - origlow
      glow=newlow
      ghigh=newhigh

!     Find the smallest of the three linewidths and use it to calculate
!     xstep, the spacing of points on the xaxis in gauss

      lwmin=LWx
      IF (LWy.lt.lwmin) lwmin=LWy
      IF (LWz.lt.lwmin) lwmin=LWz
      xstep = lwmin/10
      IF ((orighigh-origlow)/xstep .lt. 1000.)xstep=(orighigh-origlow)/1000.

      if(npts .eq. 0 .or. npts .gt. 10000)then !If we are not given the number of points we figure a useful number
        numpts=int(gps/xstep)+1
        npts = int(stickrange/xstep)-1
        IF (numpts.gt.10000) THEN
          numpts=10000
          xstep=gps/9999
        END If
      else
!        numpts = npts
        xstep = (orighigh-origlow)/(npts-1)
        numpts = int(gps/xstep)+1
      end if
      npad=int((5*lwmax)/xstep)+1
      if(numpts == 10000)npts = 10000 - 2*npad
      field(1) = origlow
      do i=2,numpts
        field(i) = field(i-1)+xstep
      end do

      ! If we are given zero values for LT and LP
      ! Estimate the spectral width and use that to determine the number of angles for the integration
      if(LT == 0 .or. LP == 0)then

        anmax = 0.
        spinmax = nucspin(1)
        do j = 1, nnuc
          Amax = Atens(j,1)
          Amin = Atens(j,1)
          do i = 1, 3
            if(Amin .gt. Atens(j,i)) Amin = Atens(j,i)
            if(Atens(j,i) .gt. Amax)then
              Amax = Atens(j,i)
!              spinmax = nucspin(j)
            end if
          end do
          if(anmax .lt. (amax-amin))then
            anmax = amax - amin
            spinmax = nucspin(j)
          end if
        end do
        gmin = g(1)
        gmax = g(1)
        do i = 1, 3
          if(gmin .gt. g(i))gmin = g(i)
          if(gmax .lt. g(i))gmax = g(i)
        end do
        specwidth = hnu/(gmin*bohrmag) - hnu/(gmax*bohrmag)
        LP = int(specwidth*2.0/lwmin + anmax*20/lwmin)
        LT = 2*LP
      end if
      if(LP .lt. 25)then
        LP = 25
        LT = 50
      end if
      allocate(res(nucmult*LP+1),istick(nucmult*LP+1))
      allocate(cop(lp),sip(lp),co2p(lp),si2p(lp),cot(lt), &
     &sit(lt),co2t(lt),si2t(lt))


!     Reset the field array
        H=0.
        stick = 0.

!     Caclulate and store the M sub I  values
      DO j=1,2
         DO i=1,nnuc
            nucqnum(i,j)=nucspin(i)
         END DO
      END DO
      DO r=2,nucmult
         DO i=nnuc,1,-1
            IF (nucqnum(i,r).gt.-nucspin(i)) THEN
                nucqnum(i,r)=nucqnum(i,r)-1
                DO j=i+1,nnuc
                   nucqnum(j,r)=nucspin(j)
                END DO
                GO TO 30
            END IF
         END DO
   30    CONTINUE
         IF (r .lt. nucmult) THEN
            DO i=1,nnuc
               nucqnum(i,r+1)=nucqnum(i,r)
            END DO
         END IF
      END DO

!     Calculate and store the theta and phi functions
      do i = 1, LP
         P = float(i)
         P = ((P-1.)/(float(LP)-1.))*pi
         cop(i) = cos(P)
         sip(i) = sin(P)
         co2p(i) = cos(2*P)
         si2p(i) = sin(2*P)
       end do

       do i = 1, LT
         P = float(i)
         P = (P-1.)/(float(LT)-1.)
         cot(i) = P
         sit(i) = sqrt(1-P**2)
         co2t = cos(2*acos(cot(i)))
         si2t = sin(2*asin(sit(i)))
       end do

!     Starting the theta and phi loops
      nres = 1
      theta: DO A=1,LT,1.
        phi: DO B=1,LP,1.
!        Calculate the direction cosines
          dcl = sit(a)*cop(b)
          dcm = sip(b)*sit(a)
          dcn = cot(a)

!        calculate the orientation dependant line width if we need it
          IF (lwiso) THEN
            LW=LWx
          ELSE
            LW=sqrt((dcl*LWx)**2+(dcm*LWy)**2+(LWz*dcn)**2)
            lpts=(10*LW)/xstep
            CALL gline (lshape, LW, xstep, lpts)
          END IF

! Calculate g effective
        geff = sqrt(dcl**2*g(1)**2 + dcm**2*g(2)**2 + dcn**2*g(3)**2)

!        Calculate the effective A values
        if(eulerflag) then
          do i = 1, nnuc
            S1 = 1/geff*(g(3)*dcn*coal(i)*sibe(i)+ &
            g(1)*dcl*(coal(i)*cobe(i)*coga(i)-sial(i)*siga(i))+ &
            g(2)*dcm*(coal(i)*cobe(i)*siga(i)+sial(i)*coga(i)))

            S2 = 1/geff*(g(3)*dcn*sial(i)*sibe(i)+ &
            g(1)*dcl*(sial(i)*cobe(i)*coga(i)+coal(i)*siga(i))- &
            g(2)*dcm*(sial(i)*cobe(i)*siga(i)-coal(i)*coga(i)))

            S3 = 1/geff*(g(3)*dcn*cobe(i)- &
            g(1)*dcl*sibe(i)*coga(i)+ &
            g(2)*dcm*sibe(i)*siga(i))
            aeff(i) = sqrt(Atens(i,1)**2*S1**2 + Atens(i,2)**2*S2**2 + Atens(i,3)**2*S3**2)
          end do
        else
          do i = 1, nnuc
             S1 = dcl*g(1)/geff
             S2 = dcm*g(2)/geff
             S3 = dcn*g(3)/geff
             aeff(i) = sqrt(Atens(i,1)**2*S1**2 + Atens(i,2)**2*S2**2 + Atens(i,3)**2*S3**2)
          end do
        end if

        newres = hnu/(geff*bohrmag)
        do r = 1, nucmult
          res(nres) = newres
          ! Compiler directives speed up the inner loop a little under ifort
          ! According to Vtune > 50% of clockticks go to 'do i = 1, nnuc' on a pentium 4
          ! ie 50% of ticks are going to branch misses on this loop
          !DEC$ NOVECTOR
          !DEC$ UNROLL(4)
          do i = 1, nnuc
            res(nres) = res(nres) + nucqnum(i,r)*aeff(i)
          end do
          nres = nres + 1
        end do
        ! Fall back to older, calculate and add line method for anisotopic lw case
        if(.not. lwiso)then
          call rad_lineadd(H,res-1,lshape,1.E0,xstep,newlow,lpts,nres-1)
          nres = 1
        end if

        END DO phi

        ! When the linewidth is isotropic we generate a "powder stick spectrum"
        if(lwiso)then
          do i = 1,nres-1
            res(i) = (res(i)-newlow)/gps
          end do
          do i = 1,nres-1
            istick(i) = res(i)*numpts
          end do
          do i = 1,nres-1
            if(istick(i) .gt. 0 .and. istick(i) .le. numpts)stick(istick(i)) = stick(istick(i)) + 1.
          end do
          nres = 1
        end if
      END DO theta

      ! Free storage
      deallocate(nucqnum,res,istick)
      deallocate(cop,sip,co2p,si2p,cot,sit,co2t,si2t)

      ! Use fourier convolution to apply a line shape to our stick spectrum
      ! The stick and convolute method is several fold faster than adding individual lines
      if(lwiso)then
        lwstick = lwx
        stick = stick/LT*LP
        call convolute(stick,numpts,lwstick,stickrange)
        H = stick
      end if

!     ***********************************************************
!     This loop normalizes the intensities to -2000 to 2000
!     ***********************************************************

      ymax=0
      ymin=0

      DO i=npad,numpts-npad
        IF (H(i).gt.ymax) ymax=H(i)
        IF (H(i).lt.ymin) ymin=H(i)
      END DO
      yave=abs(ymax-ymin)/2.
      yshift=2000.-ymax/yave*2000.

      H(1:numpts) = H(1:numpts)/yave*2000+yshift
      do i=1,numpts-npad
        H(i) = H(i+npad)
      end do
      end subroutine radsim

!     **********************************************************
!     Subroutine to add a pre-calculated line to a spectrum
!     **********************************************************

      subroutine rad_lineadd(H,HR,lshape,TP,xstep,hlow,lpts,nres)

      real :: lshape(10000), H(10000),HR(10000)

      real :: xstep,hlow,TP

      integer :: respt, midpt, step, lpts, clip,i,nres,j

      TP = 1
      midpt = lpts/2
      do j = 1, nres
        if(HR(j) .lt. hlow)cycle
        step = 0
        respt = int((HR(j)-hlow)/xstep)
        clip = 0
        if (midpt .ge. respt) clip=midpt-respt
        do i=respt-midpt+clip+1,respt+midpt
          step = step +1
          H(i) = H(i) + lshape(step+clip)
        end do
      end do
      return
      END subroutine rad_lineadd


      end module radical_simulation




