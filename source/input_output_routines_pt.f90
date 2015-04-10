module input_output
contains

!********************************************************************
! A subroutine to prompt for spectral parameters from the terminal
!********************************************************************

subroutine spec_from_terminal(specpt,pflagin,parm_prompt)
use global
use spectrum_parameters
use atom_types
use init_atoms

implicit none

! Arguments
type (spec_parms), intent(inout) :: specpt
type(prompt_flag), intent(in) :: pflagin ! A collection of logical variables to toggle prompting for values
character(len=80), optional, intent(in) :: parm_prompt

! Local variables
character(len=5) :: testnuc
character(len=50) :: nuclear_prompt
integer :: k,j,nnuc, i, limit
type (isotope) :: atom(natoms)
type(prompt_flag) :: pflag

! Make a local copy ot the prompt flag
      pflag = pflagin

!    Fill the atom list
      call fill_atoms(atom)

      if(pflag%freq) then
        write(*,*) 'Enter the microwave frequency in GHz,or -1 to quit'
        read(*,*) specpt%freq
        pflag%freq = .false.
        if(specpt%freq.lt.0)return
      end if

      if(pflag%field)then
        Write(*,*) 'Enter range for the simulation'
        read(*,*) specpt%origlow, specpt%orighigh
        pflag%field = .false.
      end if


      if(pflag%ne)then
        specpt%numelectrons = 1
        pflag%ne = .false.
        if(specpt%numelectrons == 2)then
          pflag%zfs = .true.
          pflag%intang = .true.
        else
          pflag%zfs = .false.
        end if
      end if

      if(pflag%fname) then
        write(*,*) 'Enter a filename for the simulation'
        read(*,*) specpt%fname
        pflag%fname = .false.
      endif

      if(pflag%zfs)then
        if(present(parm_prompt))write(*,'(a,f10.1,f8.1,f8.1)') trim(parm_prompt)//' ',specpt%j,specpt%d,specpt%e
        write(*,*) 'Enter J D and E'
        read(*,*) specpt%j,specpt%d,specpt%e
        write(*,*) 'Enter the inter-electron angles: theta phi psi'
        read(*,*) specpt%theta, specpt%phi, specpt%psi
        pflag%zfs = .false.
      end if

      if(pflag%intang)then
        if(present(parm_prompt))write(*,'(a,I4,I4)') trim(parm_prompt)//' ',specpt%LT, specpt%LP
        write(*,*) 'Enter the number of angles for the &
     &  integration theta phi.'
        read(*,*) specpt%LT, specpt%LP
        pflag%intang = .false.
      end if

      if(specpt%numnuclei .gt. 0)then
        limit = specpt%numnuclei
      else
        limit = 2*maxatoms
      end if

!  Loop over the 2 electrons
      nnuc = 1
      electrons: Do k=1,specpt%numelectrons
      if (pflag%gtens(k)) then
         Write (*,'(a,I1)') 'Enter gxx, gyy, gzz for electron # ',k
         if(k == 1)then
           if(present(parm_prompt))write(*,'(a,f8.4,f8.4,f8.4)') trim(parm_prompt)//' ',specpt%g1(1:3)
           Read  (*,*) specpt%g1(1:3)
         else
           if(present(parm_prompt))write(*,'(a,f8.4,f8.4,f8.4)') trim(parm_prompt)//' ',specpt%g2(1:3)
           read(*,*) specpt%g2(1:3)
         end if
         pflag%gtens(k) = .false.
      end if
      if(pflag%lw(k))then
        if(present(parm_prompt))write(*,'(a,f5.1,f5.1,f5.1)') trim(parm_prompt)//' ',specpt%lw(k,1:3)
        write (*,'(a,I1,a)') 'Enter the linewidth for electron ', k,' components(3)'
        read (*,*) specpt%lw(k,1:3)
        pflag%lw(k) = .false.
      end if

      ! This section gets information on the nuclear spin
!      nnuc = 1
      nuclear_prompt = 'Enter a nucleus, beta, equiv or done'
      nuclei: do
        if(pflag%nucspins(k))then
          write(*,*) trim(nuclear_prompt)
           read(*,*) testnuc

           ! For beta, and equiv we make multiple passes through this case construct
           select case (testnuc)
           case('done')
             pflag%nucspins(k) = .false.
             exit

           case('?')
             write(*,*) "The list of all available atom types follows"
             do i = 1, natoms
               write(*,'(I2,A)') i, ' '//trim(adjustl(atom(i)%isoname))
             end do
             nuclear_prompt = 'Enter a nucleus, beta, equiv or done'
             cycle nuclei

           case('beta')
             specpt%geom_type(nnuc) = 1
             nuclear_prompt = 'Enter the beta nucleus'
             cycle nuclei

             case default
               nuclear_prompt = 'Enter a nucleus, beta, equiv or done'

           end select

           ! If an ordinary atom name is given we come to this point
           ! Is the name on the atom list?
           do j = 1,natoms
              if(testnuc .eq. atom(j)%isoname) then
                 specpt%numnuclei = specpt%numnuclei +1
!                 nnuc = nnuc + 1
                 if(specpt%numnuclei .gt. maxatoms)exit
                 pflag%Atens(specpt%numnuclei) = .true.
                 specpt%which_elec(specpt%numnuclei) = k
                 specpt%nucname(specpt%numnuclei) = trim(atom(j)%isoname)
                 specpt%what_isotope(specpt%numnuclei)=j
                 exit
              endif
              if(j .eq. natoms)then
                write(*,*) testnuc, 'is not on my list'
                cycle nuclei
              endif
           end do
         end if

!         print*, 'numnuclei, nnuc',specpt%numnuclei, nnuc

         ! This section gets hyperfine information
         if(pflag%Atens(nnuc))then
           ! This if then elseif ladder contains the geometry appropriate prompts
           ! according to whether we have a beta, equiv, or ordinary nucleus
           if(specpt%geom_type(nnuc) == 0)then
             if(present(parm_prompt))write(*,'(a,3(f7.1),3(f6.1))') trim(parm_prompt)//' ',specpt%Atens(nnuc,1:3),&
     &        specpt%alpha(nnuc),specpt%beta(nnuc),specpt%gama(nnuc)
             write(*,'(a,a,I2)') 'Enter hyperfine: Ax Ay Az alpha beta gamma for ', &
     &        trim(adjustl(specpt%nucname(nnuc)))//' ',nnuc
             read(*,*) specpt%Atens(nnuc,1:3),specpt%alpha(nnuc),specpt%beta(nnuc), &
     &       specpt%gama(nnuc)
!             nnuc = nnuc +1
           else if(specpt%geom_type(nnuc) == 1)then
             write(*,*) 'The isotropic hyperfine A will be calculated from: A = A0 + A1*cos**2(phi)'
             write(*,*) 'Where phi is the dihedral angle between the nucleus and the pi orbital'
             write(*,*) 'A can be made anisotropic, using a correction tensor Ax, Ay, Az to be added to the isotropic value'
             if(present(parm_prompt))write(*,'(a,f6.1,5(f7.1))') trim(parm_prompt)//' ',specpt%Atens(nnuc,1:3),&
     &        specpt%alpha(nnuc),specpt%beta(nnuc),specpt%gama(nnuc)
             write(*,*) 'Enter phi, A0, A1, Ax, Ay, Az for ', &
     &        trim(adjustl(specpt%nucname(nnuc)))//' ',nnuc
             read(*,*) specpt%Atens(nnuc,1:3), specpt%alpha(nnuc),specpt%beta(nnuc), &
     &       specpt%gama(nnuc)
           else if(specpt%geom_type(nnuc) .lt. 0)then
             write(*,*) 'The isotropic hyperfine A will be calculated from: A = A0 + A1*cos**2(phi+theta)'
             write(*,*) 'Where phi is the angle between the pi orbital and ',&
             trim(adjustl(specpt%nucname(-specpt%geom_type(nnuc)))), &
             -specpt%geom_type(nnuc),' and theta is the angle between that nucleus and this one.'
             if(specpt%what_isotope(nnuc) == specpt%what_isotope(-specpt%geom_type(nnuc)))then
               if(present(parm_prompt))write(*,'(a,f6.1)') trim(parm_prompt)//' ',specpt%Atens(nnuc,1)
               write(*,*) 'Enter theta'
               read(*,*) specpt%Atens(nnuc,1)
             else
               write(*,*) 'A can be made anisotropic, using a correction tensor Ax, Ay, Az to be added to the isotropic value'
               if(present(parm_prompt))write(*,'(a,6(f7.1))') trim(parm_prompt)//' ',specpt%Atens(nnuc,1:3),&
     &          specpt%alpha(nnuc),specpt%beta(nnuc),specpt%gama(nnuc)
               write(*,*) 'Enter theta, A0, A1, Ax, Ay, Az'
               read(*,*) specpt%Atens(nnuc,1:3), specpt%alpha(nnuc),specpt%beta(nnuc), &
     &         specpt%gama(nnuc)
!               nnuc = nnuc +1
             end if
           end if

           pflag%Atens(nnuc) = .false. ! We've got the hyperfine for this nucleus so set the flag false
         end if

         nnuc = nnuc + 1

         ! Several tests to see if we are finished collecting information
         if(specpt%numnuclei == maxatoms)exit
         if(pflag%nucspins(k))cycle nuclei
         do i = 1, maxatoms
           if(specpt%which_elec(i) == k .and. pflag%Atens(i))cycle nuclei
         end do
         exit nuclei
         enddo nuclei
!        end if
      enddo electrons
      pflag%Atens = .false.
      pflag%nucspins = .false.

      return
      end subroutine spec_from_terminal

!  **********************************************************************
!  A subroutine to get a list of coupled nuclei from the terminal
!  **********************************************************************
subroutine nucspins(specpt,ielectron)

use global
use spectrum_parameters
use atom_types
use init_atoms

implicit none

! Arguments
type (spec_parms), intent(inout) :: specpt
integer, intent(in) :: ielectron

! Local variables
character(len=5) :: testnuc
character(len=50) :: nuclear_prompt, default_prompt
integer :: j, next_geom, i, limit
type(isotope) :: atom(natoms)

! Fill the atom list
call fill_atoms(atom)

default_prompt = 'Enter a nucleus, beta, equiv, ? or done'
nuclear_prompt = default_prompt
nuclei: do
   write(*,*) trim(nuclear_prompt)
   read(*,*) testnuc

   ! For beta, and equiv we make multiple passes through this case construct
   select case (testnuc)
   case('done')
      exit

   case('?')
      write(*,*) "The list of all available atom types follows"
      do i = 1, natoms
         write(*,'(I2,A)') i, ' '//trim(adjustl(atom(i)%isoname))
      end do
      nuclear_prompt = default_prompt
      next_geom = 0
      cycle nuclei

   case('beta')
      next_geom = 1
      nuclear_prompt = 'Enter the beta nucleus'
      cycle nuclei

   case('equiv')
     if(specpt%numnuclei == 0)then
        write(*,*) 'Must have at least one other nucleus defined'
        cycle nuclei
     end if
     write(*,*) 'Available nuclei'
     j = 0
     do i = 1, specpt%numnuclei
        if(specpt%geom_type(i) == 0)then
           write(*,*) i, specpt%nucname(i)
           j = i
        end if
     end do
     if(j == 0)then
        write(*,*) "Must have at least one nucleus defined before assigning it's equivalent"
        cycle nuclei
     end if
     write(*,*) 'Enter the number of the partner of this nucleus'
     read(*,*) i
     if(specpt%geom_type(i) /= 0 .or. i .gt. specpt%numnuclei .or. i .lt. 1)then
        write(*,*) 'You must enter a number from the list of nuclei above'
        nuclear_prompt = default_prompt
        cycle nuclei
     end if
     specpt%numnuclei = specpt%numnuclei +1
     if(specpt%numnuclei .gt. maxatoms)exit
     specpt%which_elec(specpt%numnuclei) = specpt%which_elec(i)
     specpt%nucname(specpt%numnuclei) = specpt%nucname(i)
     specpt%what_isotope(specpt%numnuclei)=specpt%what_isotope(i)
     specpt%geom_type(specpt%numnuclei) = i+maxatoms
     write(*,*) "Enter scale factors for the parameters of the equivalent nucleus (3 A's and 3 angles)"
     write(*,*) 'These will normally be 1.0, but can take any real value'
     read(*,*) specpt%Atens(specpt%numnuclei,1:3),specpt%alpha(specpt%numnuclei),specpt%beta(specpt%numnuclei), &
     &       specpt%gama(specpt%numnuclei)
     nuclear_prompt = default_prompt
     next_geom = 0
     cycle nuclei

!   case default
!      nuclear_prompt = default_prompt

   end select

   ! If an ordinary atom name is given we come to this point
   ! Is the name on the atom list?
   do j = 1,natoms
      if(testnuc .eq. atom(j)%isoname) then
         specpt%numnuclei = specpt%numnuclei +1
         if(specpt%numnuclei .gt. maxatoms)exit
         !pflag%Atens(specpt%numnuclei) = .true.
         specpt%geom_type(specpt%numnuclei) = next_geom
         specpt%which_elec(specpt%numnuclei) = ielectron
         specpt%nucname(specpt%numnuclei) = trim(atom(j)%isoname)
         specpt%what_isotope(specpt%numnuclei)=j
         next_geom = 0
         nuclear_prompt = default_prompt
         exit
         if(j .eq. natoms)then
            write(*,*) testnuc, 'is not on my list'
            cycle nuclei
         endif
      end if
   end do
enddo nuclei
end subroutine nucspins

! ***********************************************************************
! A subroutine to get the A tensors of a set of nuclei
! **********************************************************************
subroutine Avalues(specpt,ielectron)

use global
use spectrum_parameters
use atom_types
use init_atoms

implicit none

! Arguments
type (spec_parms), intent(inout) :: specpt
integer, intent(in) :: ielectron

! Local variables
character(len=5) :: testnuc
character(len=50) :: nuclear_prompt, default_prompt
integer :: j, nnuc, i, limit
type(isotope) :: atom(natoms)

! Fill the atom list
call fill_atoms(atom)

do nnuc = 1, specpt%numnuclei
   if(specpt%which_elec(nnuc) .ne. ielectron) cycle
   ! This if then elseif ladder contains the geometry appropriate prompts
   ! according to whether we have a beta, equiv, or ordinary nucleus
   if(specpt%geom_type(nnuc) == 0)then
      write(*,'(a,a,I2)') 'Enter hyperfine: Ax Ay Az alpha beta gamma for ', &
     &    trim(adjustl(specpt%nucname(nnuc)))//' ',nnuc
      read(*,*) specpt%Atens(nnuc,1:3),specpt%alpha(nnuc),specpt%beta(nnuc), &
     &    specpt%gama(nnuc)
   else if(specpt%geom_type(nnuc) == 1)then
      write(*,*) 'The isotropic hyperfine A will be calculated from: A = A0 + A1*cos**2(phi)'
      write(*,*) 'Where phi is the dihedral angle between the nucleus and the pi orbital'
      write(*,*) 'A can be made anisotropic, using a correction tensor Ax, Ay, Az to be added to the isotropic value'
      write(*,*) 'Enter phi, A0, A1, Ax, Ay, Az for ', &
     & trim(adjustl(specpt%nucname(nnuc)))//' ',nnuc
      read(*,*) specpt%Atens(nnuc,1:3), specpt%alpha(nnuc),specpt%beta(nnuc), &
     & specpt%gama(nnuc)
   else if(specpt%geom_type(nnuc) .lt. 0)then
      write(*,*) 'The isotropic hyperfine A will be calculated from: A = A0 + A1*cos**2(phi+theta)'
      write(*,*) 'Where phi is the angle between the pi orbital and ',&
       trim(adjustl(specpt%nucname(-specpt%geom_type(nnuc)))), &
       -specpt%geom_type(nnuc),' and theta is the angle between that nucleus and this one.'
      if(specpt%what_isotope(nnuc) == specpt%what_isotope(-specpt%geom_type(nnuc)))then
         write(*,*) 'Enter theta'
         read(*,*) specpt%Atens(nnuc,1)
      else
         write(*,*) 'A can be made anisotropic, using a correction tensor Ax, Ay, Az to be added to the isotropic value'
         write(*,*) 'Enter theta, A0, A1, Ax, Ay, Az'
         read(*,*) specpt%Atens(nnuc,1:3), specpt%alpha(nnuc),specpt%beta(nnuc), &
     &    specpt%gama(nnuc)
      end if
   end if
end do

end subroutine Avalues

!  ***********************************************************************
!  A subroutine to append a set of parameters to a log file
!  ***********************************************************************
      subroutine log_parameters(specpt,error,temperature,logfile,header)

      use global
      use spectrum_parameters

      implicit none

      type (spec_parms), pointer :: specpt
      real, optional :: error, temperature
      character*(*), optional :: logfile
      logical, optional :: header
      integer :: i
      character(len=3) :: adv
      character(len=50) :: astring

      ! Open the log file, default is "simul.log"
      if(.not. present(logfile))then
        open(unit=10,file='simul.log',access='append',status='unknown')
      else
        open(unit=10,file=trim(logfile),access='append',status='unknown')
      endif

      ! The log file is formatted as a csv (spreadsheet) file, one set of spectral parameters per line
      ! By default we write a header giving the parameter names
      adv = 'no'
      if(header .or. .not. present(header))then !The default behavior is to write a header
        if(present(error))write(10,70,advance=adv) 'error,'
        if(present(temperature))write(10,70,advance=adv) 'temp,'
        write(10,70,advance=adv) 'vnu, origlow, orighigh, '
        write(10,70,advance=adv) 'filename, '
        if(specpt%numelectrons == 2)write(10,70,advance=adv) 'J,D,E,theta,phi,psi,'
        if(specpt%numelectrons == 1 .and. specpt%numnuclei == 0)adv = 'yes'
        write(10,70,advance=adv) 'LT,LP,g1x,g1y,g1z,lw1x,lw1y,lw1z,'
        if(specpt%numnuclei == 0)adv = 'yes'
        if(specpt%numelectrons == 2)write(10,70,advance=adv) 'g2x,g2y,g2z,lw2x,lw2y,lw2z,'
        do i=1,specpt%numnuclei
          if(i == specpt%numnuclei)adv = 'yes'
          if(specpt%geom_type(i) == 0)then
            astring = ' Ax, Ay, Az, alpha, beta, gama,'
          else if(specpt%geom_type(i) == 1)then
            astring = ' phi, A0, A1, delta_Ax, delta_Ay, delta_Az,'
          else
            astring = ' theta, A0, A1, delta_Ax, delta_Ay, delta_Az,'
          end if
          write(10,77,advance=adv) 'e',specpt%which_elec(i),      &
     &        ' '//trim(specpt%nucname(i))//trim(astring)
        end do
      end if

      ! Write the parameters
      adv = 'no'
      if(present(error)) write(10,'(f12.5,a)',advance=adv) error,','
      if(present(temperature)) write(10,'(f12.5,a)',advance=adv) temperature,','
      write(10,73,advance=adv) specpt%freq, specpt%origlow, specpt%orighigh
      write(10,70,advance=adv) trim(specpt%fname)
      if(specpt%numelectrons == 2)write(10,80,advance=adv) specpt%j,specpt%d,specpt%e, &
     & specpt%theta,specpt%phi,specpt%psi
      if(specpt%numelectrons == 1 .and. specpt%numnuclei == 0)adv = 'yes'
      write(10,81,advance=adv) specpt%LT,specpt%LP,specpt%g1,specpt%lw(1,:)
      if(specpt%numnuclei == 0)adv = 'yes'
      if(specpt%numelectrons == 2)write(10,82,advance=adv) specpt%g2,specpt%lw(2,:)
      do i=1,specpt%numnuclei
         if(i == specpt%numnuclei)adv = 'yes'
         write(10,78,advance=adv) specpt%Atens(i,:), &
     &      specpt%alpha(i),specpt%beta(i),specpt%gama(i)
      end do

      close(10)
   80 format(6(',',f12.4))
   81 format(2(',',I4),6(',',f13.5))
   82 format(6(',',f13.5))
   77 format(a,I1,a)
   78 format(6(',',f13.5))
   73 format(3(f12.4,','))
   70 format(a)
      end subroutine log_parameters

!********************************************************************
! A subroutine to read a spectrum from a two column file
! If the file is a gnuplot-hostile csv it also writes a copy
! in tab separated (tsv) format
!********************************************************************

subroutine get_spectrum(specname,specx,specy,n,ierr)
  implicit none
  ! Argument variables
  character(len=*) :: specname
  real, intent(out) :: specx(10000), specy(10000)
  integer, intent(out) :: n, ierr
  ! Local variables
  real :: x, y
  character(len=80) :: inline
  integer :: comma, dot, fn_length, comma_space, i
  logical :: flag

  ! Open the file and read the data, with suitable error checking
  ierr = 0
!  print*,"specname^",trim(specname),"$"
  inquire(file=trim(specname), exist=flag)
  if(flag)then
    open(10,file=trim(specname),err=89,status='old')
    n = 0
    do
      read(10,*,end=100,err=90) x, y
      n = n+1
      specx(n) = x
      specy(n) = y
      cycle
90    write(*,*) 'error reading from file '//trim(specname)
      write(*,*) 'Possibly a problem with the file format'
      ierr = 1
100   exit
    end do
    goto 91
89  write(*,*) 'Error opening file '//trim(specname)
    ierr = 2
91  continue
    close(10)
  else
    write(*,*) "Couldn't find file: "//trim(specname)
  endif

   ! If this is a comma separated file gnuplot won't plot it properly
   ! In this case we write a new, white space separated file (.tsv)
   ! and return the name of the new file in specname
   ! A comma followed by a space will work with gnuplot though
   ! In which case we leave the file unmolested
   if(ierr == 0) then
     open(10, file=trim(specname),position='rewind')
     read(10,'(a80)') inline
     close(10)
     comma = index(inline,",")
     comma_space = index(inline,", ")
     if((comma .gt. 0 .and. comma_space == 0) .or. (comma /= comma_space))then
       fn_length = len_trim(specname)
       dot = index(specname,".",back=.true.)
!     print*, 'fn_length, dot', fn_length, dot
       if(dot .gt. 0 .and. fn_length-dot .le. 4)then ! Take a dot followed by 4 or less chars to be a removable extension
         specname = specname(1:dot)//"tsv"
       else
         specname = specname(1:fn_length)//".tsv"
       endif
       open(10,file=trim(specname))
         do i = 1, n
           write(10,*) specx(i),"	",specy(i)
         enddo
       close(10)
       write(*,*) " "
       write(*,*) 'Wrote a copy of the data into "'//trim(specname)//'", which will work better with gnuplot'
       write(*,*) " "
     endif
   endif
end subroutine get_spectrum

end module input_output
