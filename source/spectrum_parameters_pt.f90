      module spectrum_parameters
      ! This module defines the spec_parms data type
      ! and a pointer type real8_pointer that can be used to access it

      use atom_types

      ! Define a structure to hold spectral parameters, sufficient to calculate a spectrum
      Type spec_parms
         real :: freq = 9.25
         real :: origlow = 3000.
         real :: orighigh = 3500.
         character(len=80) :: fname = ' '
         real :: J = 0.
         real :: D = 0.
         real :: E = 0.
         real :: theta = 0.
         real :: phi = 0.
         real :: psi = 0.
         integer :: LT = 0
         integer :: LP = 0
         real :: g1(3) = 0.
         real :: g2(3) = 0.
         real :: lw(2,3) = 0.
         integer :: numelectrons = 0
         integer :: numnuclei = 0
         character(len=5), dimension(maxatoms) ::  nucname = ''
         real :: Atens(maxatoms,3) = 0.
         real :: alpha(maxatoms) = 0.
         real :: beta(maxatoms) = 0.
         real :: gama(maxatoms) = 0.
         integer :: which_elec(maxatoms) = 0
         integer :: what_isotope(maxatoms) = 0
         integer :: geom_type(maxatoms) = 0 ! This value changes the meaning of Atens, alpha, beta, gama
         ! Zero = tensor and eulers, 1 = beta nucleus, dihedral angle A0 A1
         ! negative = a secondary beta nucleus, greater than maxatoms = an equivalent nucleus
      End Type spec_parms

      ! This is used when we need to point to individual parameters within spec_parms
      TYPE real8_pointer
        real, POINTER :: p
      END TYPE real8_pointer

      end module spectrum_parameters


      module spectrum_parameter_handlers
      contains

      subroutine parms_to_A(spec)
      ! A subroutine to convert internal representation of beta nuclei, and equivalent nuclei
      ! which is usefull for fitting, to a more normal representation using A values and euler
      ! angles. The preferred representation for calculations.
      ! Note that this is not a reversible transformation, so fitting programs need to keep
      ! copies of the original parameters.

      use atom_types
      use init_atoms
      use spectrum_parameters
      use physical_constants

      implicit none

      type(spec_parms) :: spec, localparms
      type(isotope) :: atom(natoms)
      integer :: i, r
      real :: ang

      call fill_atoms(atom)

           ! Calculate A values for beta nuclei, and equivalent nuclei
      localparms = spec
      do i = 1, localparms%numnuclei
        select case(spec%geom_type(i))
        case(0) ! Ordinary nuclei
          cycle

        case(1) ! Beta nuclei
          ang = pi*spec%Atens(i,1)/180.
          localparms%Atens(i,1) = spec%Atens(i,2) &
     &    + spec%Atens(i,3)*cos(ang)**2 + spec%alpha(i)
          localparms%Atens(i,2) = spec%Atens(i,2) &
     &    + spec%Atens(i,3)*cos(ang)**2 + spec%beta(i)
          localparms%Atens(i,3) = spec%Atens(i,2) &
     &    + spec%Atens(i,3)*cos(ang)**2 + spec%gama(i)


        case(-maxatoms:-1) ! Secondary Beta nuclei
          ang = pi*(spec%Atens(i,1) + spec%Atens(abs(spec%geom_type(i)),1))/180.
          localparms%Atens(i,1) = spec%Atens(-spec%geom_type(i),2) &
     &    + spec%Atens(-spec%geom_type(i),3)*cos(ang)**2 &
     &    + spec%alpha(-spec%geom_type(i))
          localparms%Atens(i,2) = spec%Atens(-spec%geom_type(i),2) &
     &    + spec%Atens(-spec%geom_type(i),3)*cos(ang)**2 &
     &    + spec%beta(-spec%geom_type(i))
          localparms%Atens(i,3) = spec%Atens(-spec%geom_type(i),2) &
     &    + spec%Atens(-spec%geom_type(i),3)*cos(ang)**2 &
     &    + spec%gama(-spec%geom_type(i))
          localparms%Atens(i,1:3) = localparms%Atens(i,1:3)*atom(spec%what_isotope(i))%isognuc &
          /atom(spec%what_isotope(-spec%geom_type(i)))%isognuc

        case((maxatoms+1):(2*maxatoms+1)) ! Equivalent nuclei
          r = spec%geom_type(i) - maxatoms
          localparms%Atens(i,1:3) = spec%Atens(r,1:3) * spec%Atens(i,1:3)
          localparms%alpha(i) = spec%alpha(i) * spec%alpha(r)
          localparms%beta(i) = spec%beta(i) * spec%beta(r)
          localparms%gama(i) = spec%gama(i) * spec%gama(r)
          localparms%Atens(i,1:3) = localparms%Atens(i,1:3)*atom(spec%what_isotope(i))%isognuc &
          /atom(spec%what_isotope(spec%geom_type(i)-maxatoms))%isognuc

        end select
        if(localparms%geom_type(i) /= 0 .and. localparms%geom_type(i) .le. maxatoms)then
          localparms%alpha(i) = 0.
          localparms%beta(i) = 0.
          localparms%gama(i) = 0.
        end if
      end do
      spec = localparms

      end subroutine parms_to_A

      ! A subroutine to write simulation parameters in a semi-portable
      ! fasion.
      subroutine write_parms(parms,ierr)

      use spectrum_parameters
      use atom_types

      implicit none

      type(spec_parms) :: parms
      integer :: i,j,ierr

      ierr = 0
      open(11,file=trim(parms%fname)//'.par',status="replace", &
    &  form='unformatted',err=200)
      write(11) parms%fname
      write(11) parms%freq,parms%origlow,parms%orighigh
      write(11) parms%J,parms%D,parms%E,parms%theta,parms%phi,parms%psi
      write(11) parms%LT,parms%LP
      write(11) (parms%g1(i),i=1,3)
      write(11) (parms%g2(i),i=1,3)
      write(11) ((parms%lw(i,j),i=1,2),j=1,3)
      write(11) parms%numelectrons,parms%numnuclei
      write(11) (parms%nucname(i),i=1,maxatoms)
      write(11) ((parms%Atens(i,j),i=1,maxatoms),j=1,3)
      write(11) (parms%alpha(i),parms%beta(i),parms%gama(i),i=1,maxatoms)
      write(11) (parms%which_elec(i),i=1,maxatoms)
      write(11) (parms%what_isotope(i),i=1,maxatoms)
      write(11) (parms%geom_type(i),i=1,maxatoms)
      close(11)
      goto 210
 200  ierr = -1
 210  continue

      end subroutine write_parms

      ! A subroutine to read simulation parameters in a semi-portable
      ! fasion.
      subroutine read_parms(parms,parmfile,ierr)

      use spectrum_parameters
      use atom_types

      implicit none

      type(spec_parms) :: parms
      integer :: i,j,ierr
      character(len=*) :: parmfile

      ierr = 0
      open(11,file=trim(adjustl(parmfile)), &
    &  form='unformatted',err=100)
      read(11) parms%fname
      read(11) parms%freq,parms%origlow,parms%orighigh
      read(11) parms%J,parms%D,parms%E,parms%theta,parms%phi,parms%psi
      read(11) parms%LT,parms%LP
      read(11) (parms%g1(i),i=1,3)
      read(11) (parms%g2(i),i=1,3)
      read(11) ((parms%lw(i,j),i=1,2),j=1,3)
      read(11) parms%numelectrons,parms%numnuclei
      read(11) (parms%nucname(i),i=1,maxatoms)
      read(11) ((parms%Atens(i,j),i=1,maxatoms),j=1,3)
      read(11) (parms%alpha(i),parms%beta(i),parms%gama(i),i=1,maxatoms)
      read(11) (parms%which_elec(i),i=1,maxatoms)
      read(11) (parms%what_isotope(i),i=1,maxatoms)
      read(11) (parms%geom_type(i),i=1,maxatoms)
      close(11)
      goto 110
 100  ierr = -1
 110  continue 

      end subroutine read_parms


      end module spectrum_parameter_handlers


