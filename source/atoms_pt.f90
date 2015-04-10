      module atom_types
!********************************************************************
! These modules define a data type to hold isotopes
! along with a short list of useful values
! To add new isotopes, add them to the list at the bottom
! be sure to increse the value of natoms
! to reflect the length of the new list
!********************************************************************

! maxatoms is the maximum number of nuclear spins that a spectrum can contain
      integer, parameter :: maxatoms = 10

! natoms is the number of atoms in the list defined in module init_atoms below
      integer, parameter :: natoms = 10

! Define the isotope data type
      Type isotope
         character(len=5) :: isoname
         integer :: atomic_number
         real :: isospin
         real :: isognuc
      End Type isotope
      end module atom_types

      module init_atoms
      contains
      subroutine fill_atoms(atom)
      use atom_types
!  Define a few isotopes
      Type (isotope) :: atom(natoms)
      atom(1) = isotope('H', 1, 0.5, 5.5856948)
      atom(2) = isotope("D", 1, 1.0, 0.8574388)
      atom(3) = isotope("C13", 6, 0.5, 1.40483)
      atom(4) = isotope("N14", 7, 1.0, 0.4037637)
      atom(5) = isotope("N15", 7, 0.5, -0.5663826)
      atom(6) = isotope("O17", 8, 2.5, -0.757522)
      atom(7) = isotope("F", 9, 0.5, 5.257771)
      atom(8) = isotope("P", 15,  0.5, 2.26322)
      atom(9) = isotope("Mn", 25, 2.5, 1.3819)
      atom(10) = isotope("Co", 27, 3.5, 1.318)
      return
      end subroutine fill_atoms
      end module init_atoms

