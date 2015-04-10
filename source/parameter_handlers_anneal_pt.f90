module handlers

! Helper subroutines for the simulated annealing, and simplex
! EPR spectral fitting programs

contains

!********************************************************************
! A subroutine to get upper and lower bounds for spectral parameters
! all parameter sets are stored globally
!********************************************************************
subroutine getparms(parmflag)
   use global
   use input_output
   use spectrum_parameters
   use atom_types
   use init_atoms

   implicit none
   integer :: i, pcount, plist(15)
   character(len=20), dimension(15) :: pname
   character(len=80) :: pstring
   logical :: parmflag ! True if we have been passed a starting parameter set
   type (spec_parms), pointer :: parmpt
   type(prompt_flag) :: pflag ! logical variables to toggle prompting
   type(isotope) :: atom(natoms)

   call fill_atoms(atom)

   parmpt => parmslow

! Get some parameters that we aren't going to vary
   if(.not. parmflag)then
     write(*,*) 'Enter the microwave frequency in GHz'
     read(*,*) parmpt%freq
     write(*,*) 'Enter the number of angles for the integration, LT, LP'
     read(*,*) parmpt%LT, parmpt%LP
     write(*,*) 'How many electrons (1 or 2)'
     read(*,*) parmpt%numelectrons
   endif
   pflag%freq = .false.
   pflag%intang = .false.
   pflag%ne = .false.
   pflag%field = .false.
   pflag%fname = .false.
   if(parmpt%numelectrons == 1)pflag%zfs=.false.

   if(parmflag)then
     pflag%gtens = .false.
     pflag%lw = .false.
     pflag%nucspins = .false.
     pflag%Atens = .false.
   endif

   parmslow%fname = 'junk'
   parmshigh = parmslow
   parms = parmslow
   parmsnew = parmslow

   ! If we have a pre-defined parm set, prompt to see what to vary
   pcount = 0
   do i = 1, parmpt%numnuclei
     if(parmpt%geom_type(i) .lt. maxatoms)pcount=pcount+1
   enddo
   if(parmflag)then
     if(parmpt%numelectrons == 1)then
       pcount=2+pcount
       do i = 1, pcount-2
         if(parmpt%geom_type(i) .lt. maxatoms)pname(i) = "Atens: "//trim(adjustl(parmpt%nucname(i)))
       enddo
       pname(pcount-1)="g1"
       pname(pcount)="lw1"
     endif
     if(parmpt%numelectrons == 2)then
       pcount=5+pcount
       do i = 1, pcount-5
         if(parmpt%geom_type(i) .lt. maxatoms)pname(i) = "Atens: "//trim(adjustl(parmpt%nucname(i-5)))
       enddo
       pname(pcount-4)="g1"
       pname(pcount-3)="lw1"
       pname(pcount-2)="g2"
       pname(pcount-1)="lw2"
       pname(pcount)="zfs and angles"
    endif
    write(*,*) " "
    write(*,*) "The following table is a numbered list of parameters that can be fitted"
    do i = 1,pcount,2
      if(i+1 .le. pcount)write(*,'(I2," ",A,T30,I2," ",A)') i, trim(adjustl(pname(i))), i+1, trim(adjustl(pname(i+1)))
      if(i == pcount)write(*,'(I2," ",A)') i, trim(adjustl(pname(i)))
    enddo
    plist=pcount+1
    write(*,*) " "
    write(*,*) "Enter the numbers of the parameters you want to fit on one line followed by a slash (/)"
    read(*,*) plist
    ! Set the appropriate flags true
    if(parmpt%numelectrons == 1)then
      do i = 1,15
        if(plist(i) == pcount-1)pflag%gtens(1)=.true.
        if(plist(i) == pcount)pflag%lw(1)=.true.
        if(plist(i) .lt. pcount-1)pflag%Atens(plist(i))=.true.
      enddo
    endif
    if(parmpt%numelectrons == 2)then
      do i = 1,15
        if(plist(i) == pcount-4)pflag%gtens(1)=.true.
        if(plist(i) == pcount-3)pflag%lw(1)=.true.
        if(plist(i) == pcount-2)pflag%gtens(2)=.true.
        if(plist(i) == pcount-1)pflag%lw(2)=.true.
        if(plist(i) == pcount)pflag%zfs=.true.
        if(plist(i) .lt. pcount-4)pflag%Atens(plist(i))=.true.
      enddo
    endif
  endif

   Do
     write(*,*) ' '
     write (*,*) 'You will now be prompted for the lower limits of the spectral parameters'
     pstring="Current values:"
     if(parmflag)then
       call spec_from_terminal(parmpt,pflag,pstring)
     else
       call spec_from_terminal(parmpt,pflag)
     endif

       ! Make the prompt flag false for equivalent nuclei
       do i = 1,parmpt%numnuclei
         if(parmpt%geom_type(i) .gt. maxatoms)pflag%Atens(i) = .false.
       end do
       pflag%nucspins = .false.

     write(*,*) ' '
     write(*,*) 'Now enter the upper limits of the spectral paramters'
     write(*,*) 'Parameters that have the same upper and lower limits will not be varied'
     write(*,*) 'Triples (g, A, lw) that all have the same upper and lower limits will be varied isotropically'
     parmpt => parmshigh
     parmshigh = parmslow

     pstring="lower limit values:"
     call spec_from_terminal(parmpt,pflag,pstring)

     parmshigh%numnuclei = parmslow%numnuclei

     if(parmslow%numnuclei .ne. parmshigh%numnuclei)then
       write(*,*) 'high and low must have the same number of nuclei'
       cycle
     end if
     exit
   end do
   parms = parmslow
   parmsnew = parmslow

end subroutine getparms


!********************************************************************
! subroutine uses global upper and lower limits for epr parameters
! and figures out which ones to vary
! It returns a collection of pointers to the parameter bounds
! and to the paramters to be varied
! These pointers serve to insulate the fitting programs from the
! details of the spec_parms data structure
! This routine also prompts for information on isotopomers
! to be fit simultaineously with the main spectrum
!********************************************************************

subroutine setpointers(x,xlow,xhigh,nparm)

   use global
   use input_output
   use spectrum_parameters
   use atom_types
   use init_atoms

   implicit none
   real :: swap
   real, dimension(10000) :: isox
   type(real8_pointer), dimension(maxparms) :: x, xlow, xhigh
   integer :: nparm, i, j, ierr, inuc, l
   character(len=120) :: specname
   type(isotope) :: atom(natoms)

   call fill_atoms(atom)

! Figure out which parameters to vary, and assign pointers to them
   nparm = 0
   iso = 0
   if(parmslow%J /= parmshigh%J)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%J
     xhigh(nparm)%p => parmshigh%J
     x(nparm)%p => parms%J
     x_new(nparm)%p => parmsnew%J
   end if
   if(parmslow%D /= parmshigh%D)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%D
     xhigh(nparm)%p => parmshigh%D
     x(nparm)%p => parms%D
     x_new(nparm)%p => parmsnew%D
     nD = nparm
   end if
   if(parmslow%E /= parmshigh%E)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%E
     xhigh(nparm)%p => parmshigh%E
     x(nparm)%p => parms%E
     x_new(nparm)%p => parmsnew%E
     nE = nparm
   end if
   if(parmslow%g1(1) == parmslow%g1(2) .and. parmslow%g1(2) == parmslow%g1(3) .and. &
   parmshigh%g1(1) == parmshigh%g1(2) .and. parmshigh%g1(2) == parmshigh%g1(3) .and. &
   parmslow%g1(1) /= parmshigh%g1(1))then !vary g1 isotropically
     iso = iso+1
     niso(iso) = nparm + 1
     do i = 1,3
       nparm = nparm + 1
       xlow(nparm)%p => parmslow%g1(i)
       xhigh(nparm)%p => parmshigh%g1(i)
       x(nparm)%p => parms%g1(i)
       x_new(nparm)%p => parmsnew%g1(i)
     end do
   else
     do i = 1,3 ! Vary individual g's if needed
       if(parmslow%g1(i) /= parmshigh%g1(i))then
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%g1(i)
         xhigh(nparm)%p => parmshigh%g1(i)
         x(nparm)%p => parms%g1(i)
         x_new(nparm)%p => parmsnew%g1(i)
       end if
     end do
   end if
   if(parmslow%g2(1) == parmslow%g2(2) .and. parmslow%g2(2) == parmslow%g2(3) .and. &
   parmshigh%g2(1) == parmshigh%g2(2) .and. parmshigh%g2(2) == parmshigh%g2(3) .and. &
   parmslow%g2(1) /= parmshigh%g2(1))then !vary g2 isotropically
     iso = iso+1
     niso(iso) = nparm+1
     do i = 1,3
       nparm = nparm + 1
       xlow(nparm)%p => parmslow%g2(i)
       xhigh(nparm)%p => parmshigh%g2(i)
       x(nparm)%p => parms%g2(i)
       x_new(nparm)%p => parmsnew%g2(i)
     end do
   else
     do i = 1,3
       if(parmslow%g2(i) /= parmshigh%g2(i))then
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%g2(i)
         xhigh(nparm)%p => parmshigh%g2(i)
         x(nparm)%p => parms%g2(i)
         x_new(nparm)%p => parmsnew%g2(i)
       end if
     end do
   end if
  !Doesn't make sense to vary euler angles istropically
   if(parmslow%theta /= parmshigh%theta)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%theta
     xhigh(nparm)%p => parmshigh%theta
     x(nparm)%p => parms%theta
     x_new(nparm)%p => parmsnew%theta
   end if
   if(parmslow%phi /= parmshigh%phi)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%phi
     xhigh(nparm)%p => parmshigh%phi
     x(nparm)%p => parms%phi
     x_new(nparm)%p => parmsnew%phi
   end if
   if(parmslow%psi /= parmshigh%psi)then
     nparm = nparm + 1
     xlow(nparm)%p => parmslow%psi
     xhigh(nparm)%p => parmshigh%psi
     x(nparm)%p => parms%psi
     x_new(nparm)%p => parmsnew%psi
   end if
   do j = 1,parmslow%numelectrons
       if(parmslow%lw(j,1) == parmslow%lw(j,2) .and. parmslow%lw(j,2) == parmslow%lw(j,3) .and. &
       parmshigh%lw(j,1) == parmshigh%lw(j,2) .and. parmshigh%lw(j,2) == parmshigh%lw(j,3) .and. &
       parmslow%lw(j,1) /= parmshigh%lw(j,1))then !vary lw isotropically
         iso = iso+1
         niso(iso) = nparm+1
         do i=1,3
           nparm = nparm + 1
           xlow(nparm)%p => parmslow%lw(j,i)
           xhigh(nparm)%p => parmshigh%lw(j,i)
           x(nparm)%p => parms%lw(j,i)
           x_new(nparm)%p => parmsnew%lw(j,i)
         end do
       else
       do i=1,3
         if(parmslow%lw(j,i) /= parmshigh%lw(j,i))then
           nparm = nparm + 1
           xlow(nparm)%p => parmslow%lw(j,i)
           xhigh(nparm)%p => parmshigh%lw(j,i)
           x(nparm)%p => parms%lw(j,i)
           x_new(nparm)%p => parmsnew%lw(j,i)
         end if
     end do
     end if
   end do
   nuclei: do j = 1,parmslow%numnuclei
       ! The isotropic case
       if(parmslow%Atens(j,1) == parmslow%Atens(j,2) .and. parmslow%Atens(j,2) == parmslow%Atens(j,3) .and. &
       parmshigh%Atens(j,1) == parmshigh%Atens(j,2) .and. parmshigh%Atens(j,2) == parmshigh%Atens(j,3) .and. &
       parmslow%Atens(j,1) /= parmshigh%Atens(j,1))then !vary Atens isotropically
       iso = iso+1
       niso(iso) = nparm+1
       do i = 1,3
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%Atens(j,i)
         xhigh(nparm)%p => parmshigh%Atens(j,i)
         x(nparm)%p => parms%Atens(j,i)
         x_new(nparm)%p => parmsnew%Atens(j,i)
       end do
       cycle nuclei
       end if
       ! A secondary beta nucleus, same isotope as it's partner. eg the second of 2 beta protons
       if(parms%geom_type(j) .lt. 0)then
         if(parms%what_isotope(j) == parms%what_isotope(abs(parms%geom_type(j))))cycle nuclei
       end if
       ! The anisotropic case, also covers beta nuclei if they are solo, or the first of a group
       do i = 1,3
         if(parmslow%Atens(j,i) /= parmshigh%Atens(j,i))then
           nparm = nparm + 1
           xlow(nparm)%p => parmslow%Atens(j,i)
           xhigh(nparm)%p => parmshigh%Atens(j,i)
           x(nparm)%p => parms%Atens(j,i)
           x_new(nparm)%p => parmsnew%Atens(j,i)
         end if
       end do
      ! Euler variables may contain anisotropic components if this is a beta nucleus
      !Never makes sense to vary euler angles isotropically
       if(parmslow%alpha(j) /= parmshigh%alpha(j))then
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%alpha(j)
         xhigh(nparm)%p => parmshigh%alpha(j)
         x(nparm)%p => parms%alpha(j)
         x_new(nparm)%p => parmsnew%alpha(j)
       end if
       if(parmslow%beta(j) /= parmshigh%beta(j))then
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%beta(j)
         xhigh(nparm)%p => parmshigh%beta(j)
         x(nparm)%p => parms%beta(j)
         x_new(nparm)%p => parmsnew%beta(j)
       end if
       if(parmslow%gama(j) /= parmshigh%gama(j))then
         nparm = nparm + 1
         xlow(nparm)%p => parmslow%gama(j)
         xhigh(nparm)%p => parmshigh%gama(j)
           x(nparm)%p => parms%gama(j)
         x_new(nparm)%p => parmsnew%gama(j)
       end if
   end do nuclei


   do i=1,nparm ! Make sure high and low are in the proper order
     if(xlow(i)%p .gt. xhigh(i)%p) then
       swap = xlow(i)%p
       xlow(i)%p = xhigh(i)%p
       xhigh(i)%p = swap
     end if
   end do

   do i = 1,nparm ! Initial values set to the average of high and low
     x(i)%p = (xhigh(i)%p + xlow(i)%p)/2.
   end do

   ! Get isotopomers
   write(*,*) 'How many isotopomers? (0-5)'
   read(*,*) num_isotopomers
   do j = 1, num_isotopomers
     isotopomer(j) = parms
     prompt: do
     write(*,*) ' '
     write(*,*) 'Current list of nuclei'
     do i = 1, parms%numnuclei
       write(*,*) i, parms%nucname(i)
     end do
     write(*,*) 'Enter a number between 1 and',parms%numnuclei,' to choose a nucleus'
     read(*,*) inuc
     write(*,*) ' '
     write(*,*) 'Available isotopes'
     l = parms%what_isotope(inuc)
     do i = 1, natoms
       if(atom(l)%atomic_number == atom(i)%atomic_number) &
       write(*,*) i, atom(i)%isoname
     end do
     write(*,*) 'Enter the number of an isotope'
     read(*,*) i
     if(atom(l)%atomic_number /= atom(i)%atomic_number)then
       write(*,*) 'Not an isotope pair I recognise'
       cycle prompt
     end if
     isotopomer(j)%what_isotope(inuc) = i
     isotopomer(j)%nucname(inuc) = atom(i)%isoname
     isotopomer(j)%Atens(inuc,1:3) = parms%Atens(inuc,1:3)*(atom(i)%isognuc/atom(l)%isognuc)
     exit prompt
   end do prompt
     do
       write(*,*) ' '
       write(*,*) 'Enter the file name of a reference spectrum for this isotopomer'
       read(*,*) specname
       call get_spectrum(specname,isox,isoref(j,:),isopts(j),ierr)
       if(ierr == 1)stop
       if(ierr == 0)exit
     end do
     isotopomer(j)%origlow = isox(1)
     isotopomer(j)%orighigh = isox(isopts(j))
     write(*,*) ' '
     write(*,*) "Enter the microwave frequency for the reference spectrum"
     read(*,*) isotopomer(j)%freq
     if(isox(1) == isox(isopts(j)))stop
   end do ! isotopomers

end subroutine setpointers


!********************************************************************
! This is the "object function" called by the annealing routine
! it calculates a simulation and returns the error value
! between the simulation and the experimental spectrum
!********************************************************************
function sim_wrapper(myparms,refspec,npts,corr_err)

  use spectrum_parameters
  use input_output
  use radical_simulation
  use global
  use physical_constants
  use atom_types
  use init_atoms

  implicit none

  INTERFACE
  FUNCTION corrf (spec1, spec2, n)
  IMPLICIT none
  REAL(8) corrf (10000)
  real spec1(10000), spec2(10000)
  INTEGER n
  END FUNCTION corrf
  END INTERFACE

! Arguments
  real :: sim_wrapper
  type(spec_parms), intent(in) :: myparms
  real, intent(in) :: refspec(10000)
  integer, intent(in) :: npts
  logical, optional, intent(in) :: corr_err

! Local variables
  real, dimension(10000) :: testspec, field
  real :: myerror, gnuc_ratio, pts_ratio
  real, dimension(10000) :: autocorr, testcorr
  type(spec_parms) :: localparms
  type(isotope) :: atom(natoms)
  integer :: ispec = 1, i, j
  save ispec

  call fill_atoms(atom)

!  if(myparms%numelectrons==2)then
!    call calculate(myparms,testspec,field,npts,1)
!  else
    call radsim(myparms,testspec,field,npts)
!  end if


   ! The optional argument corr_err toggles between
   ! standard sum of squares error and
   ! error calculated as the sum of squared
   ! differences between correlation functions
   if(present(corr_err) .and. corr_err)then
     autocorr = corrf(refspec,refspec,npts)/npts
     testcorr = corrf(refspec,testspec,npts)/npts
     myerror = sum((autocorr(1:npts)-testcorr(1:npts))**2)/npts**2
   else
     myerror = sum((refspec(1:npts) - testspec(1:npts))**2)/npts**2
   end if

   ! If we are doing a simultaineous fit to isotopomers
   ! Those errors are calculated here
   do i = 1, num_isotopomers
     localparms = myparms
     localparms%freq = isotopomer(i)%freq
     ! Calculate Atensor values from the gnuc ratio
     do j = 1, localparms%numnuclei
       localparms%what_isotope(j) = isotopomer(i)%what_isotope(j)
       gnuc_ratio = atom(isotopomer(i)%what_isotope(j))%isognuc &
       /atom(localparms%what_isotope(j))%isognuc
       select case(localparms%geom_type(j))
       case(0)
         localparms%Atens(j,1:3) = localparms%Atens(j,1:3)*gnuc_ratio

        case(1)
          localparms%Atens(j,2:3) = localparms%Atens(j,2:3)*gnuc_ratio
          localparms%alpha(j) = localparms%alpha(j)*gnuc_ratio
          localparms%beta(j) = localparms%beta(j)*gnuc_ratio
          localparms%gama(j) = localparms%gama(j)*gnuc_ratio
        end select
      end do
!      if(localparms%numelectrons == 2)then
!        call calculate(localparms,testspec,field,isopts(i),1)
!      else
        call radsim(localparms,testspec,field,isopts(i))
!      end if

      pts_ratio = npts/isopts(i) ! A weighting factor if spectra have differing number of points

      ! Optional correlation based error
      if(present(corr_err) .and. corr_err)then
        autocorr = corrf(isoref(i,:),isoref(i,:),isopts(i))/isopts(i)
        testcorr = corrf(isoref(i,:),testspec,isopts(i))/isopts(i)
        myerror = myerror + pts_ratio*sum((autocorr(1:isopts(i))- &
        testcorr(1:isopts(i)))**2)/isopts(i)**2
      else
        myerror = myerror + pts_ratio*sum((isoref(i,1:isopts(i)) &
        - testspec(1:isopts(i)))**2)/isopts(i)**2
      end if
    end do

   ispec = ispec+1
   sim_wrapper = myerror

 end function sim_wrapper

end module handlers
